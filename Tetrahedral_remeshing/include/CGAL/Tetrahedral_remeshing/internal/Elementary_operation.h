// Copyright (c) 2025 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Iasonas Manolas, Jane Tournois

#ifndef CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H
#define CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <cstdlib>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_instrumentation.h>

#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <atomic>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/task_arena.h>
#endif

#include <algorithm>
#include <atomic>
#include <iterator>
#include <random>
#include <array>
#include <map>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>
#include <mutex>
#include <map>
#include <tuple>
#include <atomic>
#include <cstdio>
#include <ctime>

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#include <CGAL/Real_timer.h>
#include <cstddef>
#include <iostream>
#endif

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {


template <typename C3t3_, typename ElementType, typename ElementRange>
class Elementary_operation
{
public:
  using C3t3 = C3t3_;
  using Triangulation = typename C3t3::Triangulation;
  using Element_type = ElementType;
  using Element_range = ElementRange;

  Elementary_operation() = default;
  virtual ~Elementary_operation() = default;

  virtual Element_range get_elements(const C3t3& c3t3) const = 0;
  virtual bool execute_operation(const Element_type& e, C3t3& c3t3) = 0;
  virtual std::string operation_name() const = 0;
};

template <typename Operation>
class Elementary_operation_execution_sequential
{
public:
  using C3t3 = typename Operation::C3t3;
  using Element_range = typename Operation::Element_range;

  bool execute(Operation& op, C3t3& c3t3) const
  {
    Element_range candidates = op.get_elements(c3t3);
    if (candidates.empty())
      return false;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::size_t nb_done = 0;
    const std::size_t nb_candidates = candidates.size();
    CGAL::Real_timer timer;
    timer.start();
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
    std::size_t nb_processed = 0;
#endif
    for (const auto& element : candidates)
    {
      if (op.execute_operation(element, c3t3))
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        ++nb_done;
#endif
      }
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
      std::cout << "\r" << op.operation_name() << "... ("
                << ++nb_processed << "/" << nb_candidates << ")";
      std::cout.flush();
#endif
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    timer.stop();
    std::cout << op.operation_name() << ": " << nb_done << "/"
              << nb_candidates << " done ("
              << timer.time() << " sec)." << std::endl;
#endif
    return true;
  }
};

/**
* Picks a container type from the triangulation's concurrency tag:
* `SequentialContainer` when remeshing sequentially, `ConcurrentContainer`
* when remeshing in parallel. The sequential type is spelled out at every use
* site, so that making an operation parallelizable cannot silently change the
* container the sequential path uses.
*/
template <typename ConcurrencyTag,
          typename SequentialContainer,
          typename ConcurrentContainer>
struct Concurrency_selected_container
{
  using type = SequentialContainer;
};

#ifdef CGAL_LINKED_WITH_TBB
template <typename SequentialContainer, typename ConcurrentContainer>
struct Concurrency_selected_container<CGAL::Parallel_tag,
                                      SequentialContainer,
                                      ConcurrentContainer>
{
  using type = ConcurrentContainer;
};
#endif

template <typename ConcurrencyTag,
          typename SequentialContainer,
          typename ConcurrentContainer>
using Concurrency_selected_container_t =
  typename Concurrency_selected_container<ConcurrencyTag,
                                          SequentialContainer,
                                          ConcurrentContainer>::type;

/**
* The parallel counterpart of `Elementary_operation_execution_sequential`.
*
* An operation is parallelizable when, in addition to the `Elementary_operation`
* interface, it provides
*   - `bool lock_zone(const Element_type&, const C3t3&) const`, which locks
*     every element `execute_operation()` may touch, and
*   - `static constexpr bool requires_ordered_processing`, which says whether
*     the order of `get_elements()` carries meaning.
*
* Neither is a virtual of `Elementary_operation`: this class is a template, so
* they are found on the concrete operation. The sequential path never names
* them, and is therefore unaffected by parallelism being available.
*
* `requires_ordered_processing == true` drains a concurrent queue in the order
* `get_elements()` produced, so that threads still take the most-wanted
* elements first. `false` shuffles instead, to spread the threads over the
* triangulation and keep lock conflicts down.
*/
#ifdef CGAL_LINKED_WITH_TBB
#ifdef CGAL_TR_PHASE_WALL
// Diagnostic only: for every operation, the wall and CPU time of collecting
// its candidates and of running them, keyed by the operation's own name.
struct Op_stage_wall
{
  std::map<std::string, std::array<double,4>> by_op;   // collect w/c, run w/c
  ~Op_stage_wall()
  {
    for (const auto& [name, v] : by_op)
      std::cout << "OPSTAGE \"" << name << "\" collect=" << v[0] << "/" << v[1]
                << " run=" << v[2] << "/" << v[3] << std::endl;
  }
};
inline Op_stage_wall& op_stage_wall() { static Op_stage_wall w; return w; }
inline double op_cpu_seconds()
{
  struct timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  return double(ts.tv_sec) + 1e-9 * double(ts.tv_nsec);
}
struct Op_stage_scope
{
  std::array<double,4>& slot; int base; CGAL::Real_timer t; double c0;
  Op_stage_scope(std::array<double,4>& s, int b)
    : slot(s), base(b), c0(op_cpu_seconds()) { t.start(); }
  ~Op_stage_scope()
  { t.stop(); slot[base] += t.time(); slot[base+1] += op_cpu_seconds() - c0; }
};
#endif


#ifdef CGAL_TR_THREADTIME
// Where a worker thread's elapsed time goes, inside one parallel phase.
//
// Answers "is the busy time truly busy, or is it retrying and waiting?"
// directly, rather than through a CPU counter that cannot tell a spin from
// work. Elapsed time, not CPU time, so a thread descheduled inside a yield
// is charged for it.
struct Tr_thread_time
{
  double lock_ok = 0, apply = 0, unlock = 0, lock_fail = 0, wait = 0, worker = 0;
  unsigned long long n_ok = 0, n_fail = 0, n_wait = 0, n_did = 0;

  // `lock_zone()` split, filled in by the operations that bother to: a zone
  // is not all bookkeeping, and the interesting question is how much of it
  // is work the operation goes on to use. `z_walk` is a star gather the
  // operation consumes; `z_compute` is the operation's own computation,
  // which the smoothing has to do inside the zone because its lock is keyed
  // on the destination.
  double z_locks = 0, z_walk = 0, z_compute = 0, z_dest = 0;
  // How many zones an operation locked, against how many it turned out to
  // want. A pass that offers every vertex and decides afterwards pays the
  // whole difference.
  unsigned long long n_zone = 0, n_target = 0;
};
inline double tr_tt_now()
{
  struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
  return double(ts.tv_sec) + 1e-9 * double(ts.tv_nsec);
}
using Tr_tt_ets = tbb::enumerable_thread_specific<Tr_thread_time>;
inline Tr_tt_ets*& tr_tt_slot() { static thread_local Tr_tt_ets* p = nullptr; return p; }
inline std::atomic<Tr_tt_ets*>& tr_tt_shared()
{ static std::atomic<Tr_tt_ets*> p{nullptr}; return p; }
inline Tr_thread_time* tr_tt()
{
  Tr_tt_ets* e = tr_tt_shared().load(std::memory_order_relaxed);
  return e ? &e->local() : nullptr;
}
struct Tr_tt_report
{
  std::mutex m;
  std::vector<std::tuple<std::string, Tr_thread_time, double, int>> rows;
  ~Tr_tt_report()
  {
    std::map<std::string, std::tuple<Tr_thread_time, double, int>> agg;
    for (const auto& [name, t, wall, nt] : rows)
    {
      auto& [a, w, n] = agg[name];
      a.lock_ok += t.lock_ok; a.apply += t.apply; a.unlock += t.unlock;
      a.lock_fail += t.lock_fail; a.wait += t.wait; a.worker += t.worker;
      a.n_ok += t.n_ok; a.n_fail += t.n_fail; a.n_wait += t.n_wait;
      a.n_did += t.n_did; a.n_zone += t.n_zone; a.n_target += t.n_target;
      a.z_locks += t.z_locks; a.z_walk += t.z_walk;
      a.z_compute += t.z_compute; a.z_dest += t.z_dest;
      w += wall; n = nt;
    }
    std::fprintf(stderr,
      "%-30s %8s %8s %8s %8s %8s %8s %8s\n", "THREADTIME operation",
      "capacity", "APPLY", "lock_ok", "unlock", "LOCKFAIL", "WAIT", "unacc");
    for (const auto& [name, v] : agg)
    {
      const auto& [t, wall, nt] = v;
      const double cap = wall * nt;
      const double acc = t.apply + t.lock_ok + t.unlock + t.lock_fail + t.wait;
      std::fprintf(stderr,
        "THREADTIME %-19s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
        name.substr(0, 19).c_str(), cap, t.apply, t.lock_ok, t.unlock,
        t.lock_fail, t.wait, cap - acc);
      std::fprintf(stderr,
        "THREADTIME %-19s   %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%   ok=%llu fail=%llu waits=%llu wall=%.2f nt=%d\n",
        "", 100*t.apply/cap, 100*t.lock_ok/cap, 100*t.unlock/cap,
        100*t.lock_fail/cap, 100*t.wait/cap, 100*(cap-acc)/cap,
        t.n_ok, t.n_fail, t.n_wait, wall, nt);
      if (t.z_locks + t.z_walk + t.z_compute + t.z_dest > 0)
        std::fprintf(stderr,
          "THREADTIME %-19s   in lock_zone: spatial_locks=%.2f star_walk=%.2f compute=%.2f destination=%.2f\n",
          "", t.z_locks, t.z_walk, t.z_compute, t.z_dest);
      if (t.n_ok)
        std::fprintf(stderr,
          "THREADTIME %-19s   locked=%llu changed_something=%llu (%.1f%%)%s\n",
          "", t.n_ok, t.n_did, 100.0 * t.n_did / t.n_ok,
          t.n_zone ? "" : "");
      if (t.n_zone)
        std::fprintf(stderr,
          "THREADTIME %-19s   zones locked=%llu, wanted by the operation=%llu (%.1f%%)\n",
          "", t.n_zone, t.n_target, 100.0 * t.n_target / t.n_zone);
    }
  }
};
inline Tr_tt_report& tr_tt_report() { static Tr_tt_report r; return r; }
#endif

template <typename Operation>
class Elementary_operation_execution_parallel
{
public:
  using C3t3 = typename Operation::C3t3;
  using Element_type = typename Operation::Element_type;
  using Element_range = typename Operation::Element_range;

  bool execute(Operation& op, C3t3& c3t3) const
  {
#ifdef CGAL_TR_PHASE_WALL
    auto& slot_ = op_stage_wall().by_op[op.operation_name()];
#endif
    std::vector<Element_type> candidates;
    {
#ifdef CGAL_TR_PHASE_WALL
      Op_stage_scope scope_(slot_, 0);
#endif
      candidates = collect(op, c3t3);
    }
    if (candidates.empty())
      return false;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    CGAL::Real_timer timer;
    timer.start();
    const std::size_t nb_candidates = candidates.size();
#endif

    {
#ifdef CGAL_TR_PHASE_WALL
      Op_stage_scope scope_(slot_, 2);
#endif
#ifdef CGAL_TR_THREADTIME
      Tr_tt_ets tt_ets_;
      tr_tt_shared().store(&tt_ets_, std::memory_order_relaxed);
      const double pw0_ = tr_tt_now();
#endif
      if constexpr (Operation::requires_ordered_processing)
        run_ordered(candidates, op, c3t3);
      else
        run_unordered(candidates, op, c3t3);
#ifdef CGAL_TR_THREADTIME
      const double pwall_ = tr_tt_now() - pw0_;
      tr_tt_shared().store(nullptr, std::memory_order_relaxed);
      // The worker count comes from the threads that actually showed up, NOT
      // from tbb::this_task_arena::max_concurrency(): the enclosing arena
      // still reports the machine's width when the run was asked for fewer
      // threads, which made the capacity column -- and every percentage
      // derived from it -- four times too large on a 1-thread run.
      const int nworkers_ = static_cast<int>(tt_ets_.size());
      {
        auto& rep_ = tr_tt_report();
        std::lock_guard<std::mutex> g_(rep_.m);
        for (const Tr_thread_time& t : tt_ets_)
          rep_.rows.emplace_back(op.operation_name(), t, 0.0, 0);
        rep_.rows.emplace_back(op.operation_name(), Tr_thread_time{}, pwall_, nworkers_);
      }
#endif
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    timer.stop();
    std::cout << op.operation_name() << ": " << nb_candidates
              << " candidates (" << timer.time() << " sec, parallel)."
              << std::endl;
#endif
    return true;
  }

private:
  static std::vector<Element_type> collect(const Operation& op, const C3t3& c3t3)
  {
    Element_range range = op.get_elements(c3t3);
    if constexpr (std::is_same_v<Element_range, std::vector<Element_type>>)
      return std::move(range);
    else
    {
      std::vector<Element_type> candidates;
      candidates.reserve(std::distance(range.begin(), range.end()));
      std::copy(range.begin(), range.end(), std::back_inserter(candidates));
      return candidates;
    }
  }

  /**
  * Locks the zone of `element`, retrying until it succeeds, then runs the
  * operation and releases the zone. The retry is unbounded on purpose: a
  * failed lock means another thread holds part of this zone, and every zone
  * is released as soon as its operation ends, so the wait is finite.
  */
  static void apply_one(const Element_type& element, Operation& op, C3t3& c3t3)
  {
#ifdef CGAL_TR_LOCKCOUNT
    Lockcount_counters& lc = lockcount_counters();
    ++lc.ops;
#endif
#ifdef CGAL_TR_THREADTIME
    Tr_thread_time* tt_ = tr_tt();
    const double a_ = tt_ ? tr_tt_now() : 0.0;
#endif
    while (!op.lock_zone(element, c3t3))
    {
#ifdef CGAL_TR_LOCKCOUNT
      ++lc.retries;
#endif
      c3t3.triangulation().unlock_all_elements();
      std::this_thread::yield();
#ifdef CGAL_TR_THREADTIME
      if (tt_) ++tt_->n_wait;
#endif
    }
#ifdef CGAL_TR_THREADTIME
    const double b_ = tt_ ? tr_tt_now() : 0.0;
#endif
    op.execute_operation(element, c3t3);
#ifdef CGAL_TR_THREADTIME
    const double c_ = tt_ ? tr_tt_now() : 0.0;
#endif
    c3t3.triangulation().unlock_all_elements();
#ifdef CGAL_TR_THREADTIME
    if (tt_) { tt_->wait += b_ - a_; tt_->apply += c_ - b_;
               tt_->unlock += tr_tt_now() - c_; }
#endif
  }

  /**
  * How many elements a worker of an ORDERED operation may hold back before it
  * stops taking new ones.
  *
  * An ordered operation cannot defer to the end of the pass: its candidate
  * order carries meaning -- shortest edge first for the collapse, longest for
  * the split -- and postponing an element until everything else is done
  * changes the trajectory rather than the schedule. What it can do is stop
  * WAITING: a worker that cannot take a zone sets the element aside, takes the
  * next candidate, and comes back to the set-aside ones as soon as a few have
  * accumulated, so nothing travels more than a few places from where its
  * priority put it.
  */
  static constexpr std::size_t max_postponed = 8;

  // How many candidates a worker claims at once in run_ordered().
  //
  // The candidates are handed out in blocks rather than one at a time because
  // the hand-out point is shared by every worker: taking one element is one
  // atomic read-modify-write on a single word, and that word's cache line has
  // to travel to whichever core wants the next element. One line, one owner at
  // a time, once per element -- so the cost of handing work out grows with the
  // number of workers while the work itself does not. A block of 16 pays it
  // once per 16 elements instead.
  //
  // It stays small because the order the candidates arrive in carries meaning
  // for these operations: get_elements() puts the most-wanted first, and a
  // block is the distance by which a worker may run ahead of that order.
  static constexpr std::size_t claim_block = 16;

  static void run_ordered(std::vector<Element_type>& candidates,
                          Operation& op, C3t3& c3t3)
  {
    // A shared cursor into the candidates, not a concurrent queue: the queue
    // holds a second copy of every element and is paid per element, while the
    // cursor reads the candidates where they already are and is paid per
    // block. A worker still takes the earliest block nobody has claimed, so
    // the order elements are STARTED in is the order get_elements() produced,
    // as it was with the queue.
    //
    // At one thread the blocks are claimed 0..15, 16..31, ... by the only
    // worker, which is the candidate order exactly -- the single-threaded
    // result is unchanged, byte for byte.
    std::atomic<std::size_t> next_claim{0};
    const std::size_t n_candidates = candidates.size();

    tbb::parallel_for(0, tbb::this_task_arena::max_concurrency(),
                      [&](int)
                      {
#ifdef CGAL_TR_THREADTIME
                        Tr_thread_time* tw_ = tr_tt();
                        const double w0_ = tw_ ? tr_tt_now() : 0.0;
#endif
                        std::vector<Element_type> postponed;
                        for (;;)
                        {
                          const std::size_t first =
                            next_claim.fetch_add(claim_block, std::memory_order_relaxed);
                          if (first >= n_candidates)
                            break;
                          const std::size_t last =
                            (std::min)(first + claim_block, n_candidates);

                          for (std::size_t i = first; i < last; ++i)
                          {
                            if (!try_apply_one(candidates[i], op, c3t3))
                              postponed.push_back(candidates[i]);

                            if (postponed.size() >= max_postponed)
                              retry_postponed(postponed, op, c3t3);
                          }
                        }
                        // Whatever is still held back is taken with the
                        // waiting form: the pass is over for this worker, so
                        // there is nothing else for it to do meanwhile.
                        for (const Element_type& e : postponed)
                          apply_one(e, op, c3t3);
#ifdef CGAL_TR_THREADTIME
                        if (tw_) tw_->worker += tr_tt_now() - w0_;
#endif
                      });
  }

  // One attempt at each held-back element, in the order they were held back;
  // what still fails stays held back.
  static void retry_postponed(std::vector<Element_type>& postponed,
                              Operation& op, C3t3& c3t3)
  {
    std::size_t kept = 0;
    for (std::size_t i = 0; i < postponed.size(); ++i)
    {
      if (!try_apply_one(postponed[i], op, c3t3))
        postponed[kept++] = postponed[i];
    }
    postponed.resize(kept);
  }

  /**
  * One attempt at `element`. Unlike `apply_one()` it does not wait: a zone it
  * cannot take is given back and the element reported undone, for the caller
  * to come back to.
  */
  static bool try_apply_one(const Element_type& element, Operation& op, C3t3& c3t3)
  {
#ifdef CGAL_TR_LOCKCOUNT
    Lockcount_counters& lc = lockcount_counters();
    ++lc.ops;
#endif
#ifdef CGAL_TR_THREADTIME
    Tr_thread_time* tt_ = tr_tt();
    const double a_ = tt_ ? tr_tt_now() : 0.0;
#endif
    if(!op.lock_zone(element, c3t3))
    {
#ifdef CGAL_TR_LOCKCOUNT
      ++lc.retries;
#endif
      c3t3.triangulation().unlock_all_elements();
#ifdef CGAL_TR_THREADTIME
      if (tt_) { tt_->lock_fail += tr_tt_now() - a_; ++tt_->n_fail; }
#endif
      return false;
    }
#ifdef CGAL_TR_THREADTIME
    const double b_ = tt_ ? tr_tt_now() : 0.0;
#endif
#ifdef CGAL_TR_THREADTIME
    const bool did_ = op.execute_operation(element, c3t3);
#else
    op.execute_operation(element, c3t3);
#endif
#ifdef CGAL_TR_THREADTIME
    const double c_ = tt_ ? tr_tt_now() : 0.0;
#endif
    c3t3.triangulation().unlock_all_elements();
#ifdef CGAL_TR_THREADTIME
    if (tt_) { tt_->lock_ok += b_ - a_; tt_->apply += c_ - b_;
               tt_->unlock += tr_tt_now() - c_; ++tt_->n_ok;
               if (did_) ++tt_->n_did; }
#endif
    return true;
  }

  static std::vector<Element_type>
  gather(tbb::enumerable_thread_specific<std::vector<Element_type>>& per_thread)
  {
    std::vector<Element_type> all;
    std::size_t n = 0;
    for(const auto& v : per_thread) n += v.size();
    all.reserve(n);
    for(const auto& v : per_thread) all.insert(all.end(), v.begin(), v.end());
    return all;
  }

  /**
  * Replays the elements a pass could not take, in parallel rounds.
  *
  * The deferred elements are exactly the ones that CONFLICTED, so they tend to
  * conflict with each other, and replaying them in parallel can replay the
  * same fight: an unconditional round-based replay was measured doing four
  * times the operations on one configuration. The round therefore has to earn
  * the next one -- as soon as a round hands back more than half of what it was
  * given, the rest is done serially, where a zone cannot fail for want of
  * another thread and every element is taken exactly once.
  */
  /**
  * How much of a round has to clear for the next round to be run in parallel.
  *
  * The rule was measured at four threads, where the tail is small; at 24 it is
  * 27.6% of the wall on 94665_cdt f=0.3 (3.66 s of it in the internal flips
  * alone), so where the cutoff sits is worth a sweep rather than a constant.
  * CGAL_TR_DEFERRED_CLEARED is the fraction of `todo` a round must clear to
  * earn another parallel one, in percent. 50 is the rule as measured at four
  * threads and is the default; 0 lets any progress at all earn another round,
  * which is the fully parallel replay; 100 can only be met by a round that
  * clears everything, so the leftovers always go to the serial form. Read once
  * per process, never on a locking path.
  */
  static int deferred_cleared_percent()
  {
    static const int pct = []
    {
      if(const char* const env = std::getenv("CGAL_TR_DEFERRED_CLEARED"))
      {
        const int n = std::atoi(env);
        if(n >= 0 && n <= 100)
          return n;
      }
      return 50;
    }();
    return pct;
  }

  static void run_deferred(std::vector<Element_type> todo,
                           Operation& op, C3t3& c3t3)
  {
    const int cleared_pct = deferred_cleared_percent();
    while(!todo.empty())
    {
      tbb::enumerable_thread_specific<std::vector<Element_type>> again;
      tbb::parallel_for_each(todo,
                             [&](const Element_type& element)
                             {
                               if(!try_apply_one(element, op, c3t3))
                                 again.local().push_back(element);
                             });
      std::vector<Element_type> next = gather(again);
      if(next.empty())
        return;
      // The round cleared `todo.size() - next.size()`; it earns another one
      // only if that is at least `cleared_pct` percent of what it was given.
      if(100 * (todo.size() - next.size()) < std::size_t(cleared_pct) * todo.size())
      {
        // ONE thread finishes what the rounds could not, with the waiting
        // form of apply_one. This is a serial section, it is not small, and it
        // GROWS with the thread count -- measured at 4 threads on
        // 409635_cdt_0.5, one thread finishes 88 015 elements in 0.42 s of an
        // 11.3 s run, 3.7% of the whole, against nothing at all at one thread,
        // where no element is ever deferred.
        //
        // Spreading it over the threads was tried and does NOT pay: as
        // tbb::parallel_for_each over the same waiting form it measured
        // +0.572% wall (5 balanced blocks, sd 1.102%) on 1146193_cdt_1.5 at 4
        // threads, against a predicted -3.7%. These are precisely the elements
        // the rounds could not clear, so they conflict with each other; moved
        // onto several threads they spend the saving waiting in
        // std::this_thread::yield() instead. Whether that still holds at 24
        // threads, where the tail is larger, is not measurable on a 4-core
        // box -- CGAL_TR_TOPSTAGE reports its size and its wall time so that
        // it can be read off a bigger machine.
#ifdef CGAL_TR_TOPSTAGE
        CGAL_TR_TOPSTAGE_SCOPE((std::string("~serial deferred tail: ")
                                + op.operation_name()).c_str());
        top_stage_times().by_stage[std::string("~serial deferred tail elements: ")
                                   + op.operation_name()][2] += double(next.size());
#endif
        for(const Element_type& element : next)
          apply_one(element, op, c3t3);
        return;
      }
      todo.swap(next);
    }
  }

  static void run_unordered(std::vector<Element_type>& candidates,
                            Operation& op, C3t3& c3t3)
  {
    // No shuffle. It was introduced to spread threads over the mesh, and it
    // does not pay for itself: removing it is -2.807% wall time on the frozen
    // Tier-A 24 (CLEAR, 20/24 configs faster, instructions flat at +0.066%),
    // with the gain concentrated on the heavy meshes -- 1146193_cdt_1.5
    // -15.40%, 65617_cdt_0.5 -7.59%. Same work, better order: get_elements()
    // already produces candidates in an order the shuffle was destroying.
    //
    // It was also the largest single source of variance in the measurement
    // rig. Reseeding from random_device on every call makes every run remesh a
    // different sequence, which no replication inside one screen can average
    // out. A/A on the same binary, per-config wall sd falls 7.192% -> 1.401%
    // and the instruction null +0.827% -> -0.043%. One thread becomes
    // deterministic, which is what made an operation-level change measurable
    // at all.
    //
    // A worker that cannot take an element's zone does NOT wait for it. It
    // sets the element aside and takes the next one; what is set aside is
    // replayed once the pass has joined. Waiting was measured as the larger
    // half of the parallel path's cost -- 62-74% of the extra instructions a
    // 4-thread run executes over a 1-thread run are kernel instructions, and
    // they are `sched_yield()` in the retry loop, which almost never finds
    // another runnable task to switch to.
    //
    // Only the UNORDERED operations may do this. An ordered one takes its
    // candidate order from get_elements() for a reason, and deferring there
    // changes the trajectory: deferring every operation was measured running
    // 6.6 M operations against 4.55 M at one thread.
    //
    // At one thread nothing is ever deferred -- `lock_zone()` cannot fail when
    // no other thread holds anything -- so the single-threaded result is
    // exactly what it was.
    tbb::enumerable_thread_specific<std::vector<Element_type>> deferred;
    tbb::parallel_for_each(candidates,
                           [&](const Element_type& element)
                           {
                             if(!try_apply_one(element, op, c3t3))
                               deferred.local().push_back(element);
                           });
    run_deferred(gather(deferred), op, c3t3);
  }
};
#endif // CGAL_LINKED_WITH_TBB

/**
* Selects the execution strategy from the triangulation's concurrency tag.
* `Elementary_operation_executor<Op, Tag>` is the executor to instantiate.
*/
template <typename Operation, typename ConcurrencyTag>
struct Elementary_operation_executor_selector
{
  using type = Elementary_operation_execution_sequential<Operation>;
};

#ifdef CGAL_LINKED_WITH_TBB
template <typename Operation>
struct Elementary_operation_executor_selector<Operation, CGAL::Parallel_tag>
{
  using type = Elementary_operation_execution_parallel<Operation>;
};
#endif

template <typename Operation, typename ConcurrencyTag>
using Elementary_operation_executor =
  typename Elementary_operation_executor_selector<Operation, ConcurrencyTag>::type;

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H
