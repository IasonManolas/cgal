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

// for next_wave_claim_stamp(), used by the conflict-free wave path below.
// Required for GCC/Clang two-phase lookup : it is a non-dependent name, so it
// must be declared before the template that calls it.
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <string>

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#include <CGAL/Real_timer.h>
#include <cstddef>
#include <iostream>
#endif

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/task_group.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_sort.h>
#include <tbb/combinable.h>
#include <tbb/concurrent_vector.h>
#include <tbb/task_arena.h>
#include <functional>
#include <cstdlib>
#include <atomic>
#include <iostream>
#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/functional/hash.hpp>
#endif

#ifndef LOCK_GRID_SIZE
#  define LOCK_GRID_SIZE 128 // grid cells per axis for the spatial lock grid
#endif

// Sub-phase wall-clock breakdown of the unordered-processing pipeline
// (candidate collection, kd-partition, the classification passes, the apply
// pass). Diagnostic only: prints one line per operation per call.
#ifdef CGAL_TR_SUBPHASE_TIMINGS
#  include <chrono>
#  include <iostream>
#  define CGAL_TR_SUB_VARS()      std::chrono::steady_clock::time_point sub_t0_;      double sub_[8] = {0., 0., 0., 0., 0., 0., 0., 0.}
#  define CGAL_TR_SUB_T0() sub_t0_ = std::chrono::steady_clock::now()
#  define CGAL_TR_SUB_ADD(i)      sub_[i] += std::chrono::duration<double>(        std::chrono::steady_clock::now() - sub_t0_).count()
#else
#  define CGAL_TR_SUB_VARS() ((void)0)
#  define CGAL_TR_SUB_T0()   ((void)0)
#  define CGAL_TR_SUB_ADD(i) ((void)0)
#endif

#ifdef CGAL_TETRAHEDRAL_REMESHING_LOCK_STATS
// How many unordered-processing elements took the lock-free interior path, how
// many had to lock, and how many lock attempts failed. Diagnostic only : the
// counters are relaxed atomics, so they perturb the very contention they
// measure as little as possible, but they are still not free -- do not enable
// them in a timing run.
namespace lock_stats {
  inline std::atomic<std::size_t>& interior()
  { static std::atomic<std::size_t> c{0}; return c; }
  inline std::atomic<std::size_t>& locked()
  { static std::atomic<std::size_t> c{0}; return c; }
  inline std::atomic<std::size_t>& lock_retry()
  { static std::atomic<std::size_t> c{0}; return c; }
  inline void dump()
  {
    const std::size_t i = interior().load(), l = locked().load(), r = lock_retry().load();
    const std::size_t t = i + l;
    std::cout << "[locks] interior " << i << " (" << (t ? 100.0 * i / t : 0.0)
              << "%), locked " << l << ", retries " << r << std::endl;
  }
}
#  define CGAL_TR_COUNT(which) \
     lock_stats::which().fetch_add(1, std::memory_order_relaxed)
#else
#  define CGAL_TR_COUNT(which) ((void)0)
#endif

#ifndef LB_BUCKETS_PER_THREAD
// Unordered ops partition the candidate set into ~LB_BUCKETS_PER_THREAD * nthreads
// equal-count, spatially-compact kd-buckets so TBB work-stealing balances load.
//
// One bucket per thread, not several. More buckets balance load better, but
// every extra bucket adds interface area, and only elements within one cell-ring
// of an interface have to lock : the rest run lock-free. That fraction is what
// dominates. Measured on bear.mesh at 8 threads, the elements taking the
// lock-free interior path go 25.3% at 4 buckets/thread -> 34.4% at 2 -> 56.3%
// at 1, and flip and smooth track it (flip 20.98 s -> 20.14 -> 19.46), while
// collapse, which this knob cannot affect, stays flat. Lock *contention* is not
// the issue -- only 15.6% of lock attempts ever retry -- the cost is acquiring
// 2.2 million locks at all.
//
// Overridable at run time via CGAL_TR_BUCKETS_PER_THREAD; raise it if a mesh
// turns out to be imbalanced enough that idle threads cost more than the locks.
#  define LB_BUCKETS_PER_THREAD 1
#endif

#ifndef COLLAPSE_WAVE_FACTOR
// Collapse selects conflict-free waves (see collapse_short_edges.h): only a
// fraction of the candidates scanned can join one, so its waves are sized
// well above split's to keep every thread busy.
#  define COLLAPSE_WAVE_FACTOR 256
#endif

#ifndef ORDERED_WAVE_FACTOR
// Ordered ops (split/collapse) process the globally-sorted candidate list in
// consecutive waves of ORDERED_WAVE_FACTOR * nthreads elements: parallel
// within a wave, strictly sequential across waves, so the greedy global
// priority order (e.g. longest-first) is preserved between waves.
#  define ORDERED_WAVE_FACTOR 8
#endif

#include <atomic>
#include <algorithm>
#include <vector>
#include <boost/container/small_vector.hpp>
#include <thread>
#include <chrono>
#include <utility>
#include <fstream>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

#ifdef CGAL_TR_FLIP_LOCK_PROBE
// [DEBUG-lk01] Which protection regime the current thread is executing an
// element under. Read by the lock-ownership probe in flip_edges.h so an
// unlocked read can be attributed to the pass that permitted it:
//   1 = INTERIOR: bucket lock elision, no locks held BY DESIGN. An unlocked
//       read here is only a bug if the bucket's disjointness guarantee does
//       not actually cover the footprint the read touches.
//   2 = LOCKED: lock_zone_precomputed() succeeded, so every vertex the
//       operation may touch is supposed to be held. Any unlocked read here
//       is unambiguously a hole in the lock set.
//   3 = WAVE: conflict-free wave claim.
inline int& tr_probe_mode() { static thread_local int m = 0; return m; }
struct Tr_probe_mode_scope
{
  int saved;
  explicit Tr_probe_mode_scope(int m) : saved(tr_probe_mode()) { tr_probe_mode() = m; }
  ~Tr_probe_mode_scope() { tr_probe_mode() = saved; }
};
#define CGAL_TR_PROBE_MODE(m) \
  ::CGAL::Tetrahedral_remeshing::internal::Tr_probe_mode_scope probe_mode_scope_(m)
#else
#define CGAL_TR_PROBE_MODE(m) do {} while(0)
#endif

template <typename C3t3_, typename ElementType, typename ElementRange>
class Elementary_operation
{
public:
  using C3t3 = C3t3_;
  using Triangulation = typename C3t3::Triangulation;
  using Element_type = ElementType;
  using Element_range = ElementRange;
  using Vertex_handle = typename Triangulation::Vertex_handle;

  Elementary_operation() = default;
  virtual ~Elementary_operation() = default;

  virtual Element_range get_elements(const C3t3& c3t3) const = 0;
  virtual bool execute_operation(const Element_type& e, C3t3& c3t3) = 0;

  // Bucket-scoped locking (opt-in; see flip_edges.h).
  //
  // Operations that acquire vertex locks on demand during execute_operation()
  // cannot cache anything derived from the triangulation across an unlock:
  // once the lock is gone another worker may reshape that neighbourhood. With
  // per-element unlocking that forces a cache flush per element, which is
  // expensive for flip (~18% on bear/1-thread) because consecutive elements in
  // a bucket overlap heavily.
  //
  // A bucket is processed start-to-finish by a single worker, so holding the
  // locks for the whole bucket keeps every cached entry valid for as long as
  // it is useful, and unlocking once per bucket instead of once per element.
  // This cannot deadlock: acquisition is always try-lock, and a failure
  // abandons the element (and releases everything) rather than waiting.
  virtual bool prefers_bucket_scoped_locks() const { return false; }

  // Called when the worker releases the locks it held for a bucket, so the
  // operation can drop anything cached under them.
  virtual void on_locks_released() {}
  virtual std::string operation_name() const = 0;

  // Hooks used only by Elementary_operation_execution_parallel. Defaulted so
  // existing operations and sequential builds are unaffected.
  // True if the operation's footprint extends beyond the 1-rings of its
  // locked_vertices(), so the lock-elision classification needs the extra
  // one-cell-layer halo (see apply_unordered_processing, Pass 2b). Flip walks and
  // rewrites the neighbourhoods of the opposite/ring vertices around its edge ->
  // true. Smooth only touches its own vertex's 1-ring, which the endpoint tag
  // already covers -> false, so it keeps eliding on the un-expanded tag.
  virtual bool footprint_needs_halo() const { return false; }

  virtual bool requires_ordered_processing() const { return false; }
  virtual bool lock_zone(const Element_type&, const C3t3&) const { return true; }

  // Boundary path of apply_unordered_processing's Pass 3 only: same contract
  // as lock_zone() (try-lock everything needed, return false and let the
  // caller unlock-all-and-retry on partial failure), but given the exact
  // vertex set locked_vertices() already computed for this element back in
  // Pass 1 -- before any parallel mutation started, so that computation was
  // safe to do by walking the triangulation with no lock held. An operation
  // whose true footprint needs more than {the element's own endpoints} (e.g.
  // Internal_edge_flip_operation's ring vertices) can override this to lock
  // every vertex in `lv` instead of rediscovering them here: rediscovering a
  // variable-sized vertex set live, inside Pass 3, is unsafe (it requires
  // reading the triangulation around the element while other buckets are
  // concurrently mutating unrelated parts of it) -- confirmed by a SIGSEGV
  // inside Triangulation_data_structure_3::is_edge() when an earlier attempt
  // tried exactly that. Using the Pass-1-computed `lv` sidesteps the problem
  // by never re-reading the triangulation to find out what to lock; Pass 3
  // only ever *locks* handles it already safely knows about. Default:
  // delegates to lock_zone(), unchanged for operations that have not been
  // widened this way.
  virtual bool lock_zone_precomputed(const Element_type& e, const C3t3& c3t3,
                                      const boost::container::small_vector<Vertex_handle, 2>&) const
  {
    return lock_zone(e, c3t3);
  }

  // True to route this operation through apply_conflict_free_waves() (see its
  // long comment below) rather than apply_unordered_processing()'s
  // bucket/lock-elision scheme, regardless of the CGAL_TR_UNORDERED_WAVES env
  // var. Per-operation rather than global: apply_conflict_free_waves()'s claim
  // is exactly as strong as locked_vertices() makes it (it stamps the full
  // star of every vertex locked_vertices() returns), computed serially during
  // wave selection before any parallel mutation starts -- so an operation
  // whose locked_vertices() has been widened to cover its true footprint (see
  // Internal_edge_flip_operation) gets a claim that is sound by construction,
  // with no analogue of lock_zone()'s "must already hold a lock to safely
  // read what needs locking" bootstrapping problem. Default false so
  // operations that have not had their locked_vertices() audited this way
  // keep the pre-existing bucket-based path.
  virtual bool prefers_conflict_free_waves() const { return false; }

  // True to defer this operation's Pass-3 BOUNDARY elements (see
  // apply_unordered_processing) to a conflict-free-wave selection restricted
  // to just that subset -- run once, after all of this phase's INTERIOR
  // elements have finished executing -- instead of processing each boundary
  // element individually through lock_zone()/lock_zone_precomputed(). Unlike
  // prefers_conflict_free_waves() (which routes the WHOLE candidate set
  // through a serial wave selection, replacing apply_unordered_processing()
  // entirely -- correct but too slow for a dense candidate graph, since claim
  // overlap is common across the whole mesh), this keeps the cheap
  // interior/boundary split and only pays the serial wave-selection cost for
  // the boundary subset, which is a small minority of candidates by
  // construction. Sound for the same reason apply_conflict_free_waves() is:
  // the wave claim stamps the full star of every locked_vertices() vertex,
  // computed from the Pass-1-cached `all_lv` (no re-discovery), and boundary
  // elements only ever start running after every interior element in this
  // phase has finished -- interior elements by definition never touch a
  // boundary-tagged vertex, so nothing about their execution can have
  // invalidated the triangulation state the boundary wave selection reads.
  // Default false so operations that have not been audited this way keep the
  // pre-existing per-element lock_zone() path.
  virtual bool prefers_boundary_waves() const { return false; }

  // For prefers_boundary_waves() + footprint_needs_halo() operations only:
  // how many of locked_vertices()'s LEADING entries to skip when growing the
  // wave claim's extra halo layer (see apply_unordered_processing's boundary
  // wave pass). Default 0: grow from every vertex in lv. Override to the
  // number of "primary" vertices at the front of lv (e.g. 2 for an edge's own
  // {e.first, e.second}) when those specific vertices' own incident-cell
  // stars are not local to the element's true footprint -- growing from them
  // pulls in every cell touching that vertex anywhere in the mesh, not just
  // the cells near this element, which showed up as a severe over-claim
  // (spurious conflicts, effectively a livelock) for Internal_edge_flip_operation.
  virtual std::size_t halo_growth_skip() const { return 0; }

  // The vertices whose incident-cell 1-rings lock_zone() locks. Used by
  // apply_unordered_processing for interior/boundary lock elision: an element all of
  // whose locked vertices are bucket-interior (never shared with another
  // concurrently-processed bucket) can run lock-free. Default: empty -> the element
  // is conservatively treated as BOUNDARY (always locked). Flip and smooth override.
  virtual void locked_vertices(const Element_type&, const C3t3&,
                                boost::container::small_vector<Vertex_handle, 2>&) const {}

  // Buckets per thread for apply_unordered_processing's kd-partition, for
  // operations whose optimum differs from the global default. 0 (the default)
  // means "use the global LB_BUCKETS_PER_THREAD".
  //
  // The knob trades two costs that do NOT move together across operations.
  // More buckets add interface area, and only elements near an interface have
  // to lock -- that is what LB_BUCKETS_PER_THREAD=1 was tuned against, and it
  // is still right for smooth, which is bound by the number of locks acquired.
  // Flip is not: it holds bucket-scoped locks, so an interface element cannot
  // acquire until the *neighbouring bucket finishes*, and with one bucket per
  // thread that is a ~100k-element wait. Flip therefore wants many small
  // buckets and smooth wants few, and a single global value cannot serve both.
  // See the flip bucket sweep in the parallel optimization log.
  virtual std::size_t buckets_per_thread_hint() const { return 0; }
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

#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB

template <typename Operation>
class Elementary_operation_execution_parallel
{
public:
  using C3t3 = typename Operation::C3t3;
  using Element_type = typename Operation::Element_type;
  using Element_range = typename Operation::Element_range;
  using Vertex_handle = typename C3t3::Vertex_handle;

  bool execute(Operation& op, C3t3& c3t3) const
  {
    CGAL_TR_SUB_VARS();
    CGAL_TR_SUB_T0();
    std::vector<Element_type> candidates = collect_candidates(op, c3t3);
    CGAL_TR_SUB_ADD(0);
#ifdef CGAL_TR_SUBPHASE_TIMINGS
    std::cout << "[sub] " << op.operation_name()
              << " | get_elements " << sub_[0]
              << " s | candidates " << candidates.size() << std::endl;
#endif
    if (candidates.empty())
      return false;

    ensure_lock_data_structure_initialized(c3t3);

    if (op.requires_ordered_processing())
      return apply_ordered_processing(candidates, op, c3t3);
    else if (op.prefers_conflict_free_waves() || use_unordered_waves())
      return apply_conflict_free_waves(candidates, op, c3t3);
    else
      return apply_unordered_processing(candidates, op, c3t3);
  }

private:
  // Runtime rather than compile-time switch so the two unordered strategies can
  // be measured alternately within a single build. Timings on this workload
  // drift by ~20% over an hour of sustained load, which is larger than the
  // effect being measured, so A and B have to be interleaved rather than run as
  // separate batches -- and that is impossible if choosing between them needs a
  // rebuild. Read once; the default is the kd-bucket path.
  // Overridable at run time for the same reason as use_unordered_waves(): this
  // trades load balance against lock cost -- more buckets balance better but
  // add bucket-interface area, and only interface elements pay for lock_zone --
  // so the optimum has to be found by interleaved A/B, not batched rebuilds.
  // 0 when the environment does not set it, so the per-operation hint can be
  // consulted. An explicit CGAL_TR_BUCKETS_PER_THREAD still wins over the hint,
  // which is what makes a global sweep across all operations possible.
  static std::size_t buckets_per_thread_env()
  {
    static const std::size_t n = []() -> std::size_t
      {
        const char* const e = std::getenv("CGAL_TR_BUCKETS_PER_THREAD");
        if (e != nullptr)
        {
          const long v = std::strtol(e, nullptr, 10);
          if (v > 0)
            return static_cast<std::size_t>(v);
        }
        return std::size_t(0);
      }();
    return n;
  }

  static std::size_t buckets_per_thread(const Operation& op)
  {
    if (const std::size_t e = buckets_per_thread_env())
      return e;
    if (const std::size_t h = op.buckets_per_thread_hint())
      return h;
    return static_cast<std::size_t>(LB_BUCKETS_PER_THREAD);
  }

  static bool use_unordered_waves()
  {
    static const bool on = []
      {
        const char* const e = std::getenv("CGAL_TR_UNORDERED_WAVES");
        return e != nullptr && *e != '0';
      }();
    return on;
  }

  // Unordered processing, conflict-free variant. Instead of letting threads
  // collide and arbitrating with lock_zone(), build each wave so its members
  // provably cannot interfere, then run the wave with no locking at all.
  //
  // Selecting an element stamps every vertex of every cell incident to its
  // locked_vertices(); a later candidate joins the wave only if none of its
  // locked_vertices() already carries this wave's stamp. That test is sound
  // because a cell shared between a candidate and an already-selected element
  // would contain one of the candidate's locked vertices and be incident to the
  // selected element, so that vertex would already be stamped -- testing the
  // handful of locked vertices is therefore exactly as strong as walking the
  // full star, at O(1) instead of O(star), and confines the star walk to
  // accepted elements.
  //
  // Rejected candidates are carried to the next wave rather than retried in
  // place, so no thread ever spins. Waves are not size-capped: unlike the
  // ordered path there is no priority order to preserve, so a wave may take
  // every candidate that fits, which keeps the number of (serial) selection
  // passes down -- this is greedy graph colouring, and the pass count behaves
  // like a chromatic number, not like the candidate count.
  //
  // The motivation is that lock_zone() here is memory-bound, not contended:
  // profiling attributes its cost to the first touch of each vertex's point and
  // to the per-thread lock bitmap, both cache misses, paid for every candidate
  // whether or not it ends up modifying anything.
  //
  // SOUNDNESS CONDITION -- off by default, and not usable by every operation.
  // What the selection above guarantees is that no two members of a wave share
  // a *cell*. That is enough only for an operation whose reads and writes stay
  // within the cells of its own region. The flip operations do not qualify:
  // they share a mutable incident-cell cache and clear its entries for all four
  // vertices of every cell they touch (flip_edges.h), so two members with
  // disjoint cells but one vertex in common can still collide -- one clears an
  // entry while the other iterates it. Making the claim cover that would mean
  // stamping the 2-ring of each element rather than the 1-ring (that rule is
  // sound: an element sharing a vertex with region(A) must have an endpoint
  // within two cell-steps of A's, hence stamped), but the 2-ring walk costs
  // roughly fifteen times the 1-ring one, serially, per accepted element --
  // more than the locking it would replace. Hence the default below.
  bool apply_conflict_free_waves(std::vector<Element_type>& elements,
                                 Operation& op, C3t3& c3t3) const
  {
    using Cell_handle = typename C3t3::Triangulation::Cell_handle;
    auto& tr = c3t3.triangulation();

    std::vector<Element_type> pending(std::make_move_iterator(elements.begin()),
                                      std::make_move_iterator(elements.end()));
    std::vector<Element_type> wave, deferred;
    boost::container::small_vector<Vertex_handle, 2> lv;
    boost::container::small_vector<Cell_handle, 64> inc_cells;

    while (!pending.empty())
    {
      const std::size_t wave_stamp = next_wave_claim_stamp();
      wave.clear();
      deferred.clear();

      for (const Element_type& element : pending)
      {
        lv.clear();
        op.locked_vertices(element, c3t3, lv);

        // No claim information: fall back to running it alone, after the
        // parallel waves, rather than guessing at its conflict domain.
        if (lv.empty())
        {
          deferred.push_back(element);
          continue;
        }

        bool conflicts = false;
        for (const Vertex_handle& v : lv)
          if (v->wave_claim_stamp() == wave_stamp) { conflicts = true; break; }
        if (conflicts)
        {
          deferred.push_back(element);
          continue;
        }

        for (const Vertex_handle& v : lv)
        {
          inc_cells.clear();
          tr.incident_cells(v, std::back_inserter(inc_cells));
          for (const Cell_handle& c : inc_cells)
            for (int k = 0; k < 4; ++k)
              c->vertex(k)->set_wave_claim_stamp(wave_stamp);
        }
        wave.push_back(element);
      }

      // Nothing could be claimed (every element had an empty conflict domain):
      // run the remainder sequentially rather than loop forever.
      if (wave.empty())
      {
        for (const Element_type& element : deferred)
          op.execute_operation(element, c3t3);
        return true;
      }

      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, wave.size()),
        [&](const tbb::blocked_range<std::size_t>& r)
        {
          for (std::size_t i = r.begin(); i != r.end(); ++i)
            op.execute_operation(wave[i], c3t3);
        });

      pending.swap(deferred);
    }
    return true;
  }

  std::vector<Element_type> collect_candidates(const Operation& op, const C3t3& c3t3) const
  {
    Element_range elements = op.get_elements(c3t3);
    std::vector<Element_type> candidates;
    if constexpr (std::is_same<Element_range, std::vector<Element_type>>::value)
    {
      candidates = std::move(elements); // directly move if vector
    }
    else
    {
      tbb::combinable<std::vector<Element_type>> local_candidates;
      tbb::parallel_for_each(elements.begin(), elements.end(),
        [&](const Element_type& element) { local_candidates.local().push_back(element); });
      local_candidates.combine_each([&](const std::vector<Element_type>& local)
        { candidates.insert(candidates.end(), local.begin(), local.end()); });
    }
    return candidates;
  }

  // Ensure the lock data structure is initialized when needed
  void ensure_lock_data_structure_initialized(C3t3& c3t3) const
  {
    auto& triangulation = c3t3.triangulation();
    if (!triangulation.get_lock_data_structure())
    {
      static typename C3t3::Triangulation::Lock_data_structure lock_ds(
        c3t3.bbox(), LOCK_GRID_SIZE);
      triangulation.set_lock_data_structure(&lock_ds);
    }
  }

  // Ordered processing (Edge_split, Edge_collapse): `elements` arrives already
  // sorted in the operation's required global priority order (e.g.
  // longest-first for split). Partitioning that sorted range into
  // spatially-compact buckets and running all buckets concurrently (as
  // apply_unordered_processing does) discards the global order: an element
  // near the front of the sort can land in a bucket that a busy thread only
  // reaches after another thread has already processed a bucket containing
  // elements far later in the sort. For split/collapse this materially
  // changes which elements get processed and in what mesh state, degrading
  // quality relative to a strictly sequential run.
  //
  // Instead, process the sorted range in consecutive waves of
  // ORDERED_WAVE_FACTOR * nthreads elements: within a wave, elements run
  // concurrently (conflicts arbitrated by lock_zone as before); waves
  // themselves run strictly one after another. This keeps the greedy global
  // priority order intact across waves (the only place it matters -- wave
  // boundaries are far apart in the mesh already, so intra-wave scheduling
  // noise is confined to a small, geometrically-arbitrary neighbourhood of
  // the sort) while still parallelizing the bulk of the work.
  bool apply_ordered_processing(std::vector<Element_type>& elements, Operation& op, C3t3& c3t3) const
  {
    const std::size_t nthreads =
      static_cast<std::size_t>((std::max)(1, tbb::this_task_arena::max_concurrency()));
    const std::size_t wave_size =
      (std::max)(std::size_t(1), std::size_t(ORDERED_WAVE_FACTOR) * nthreads);

    for (std::size_t lo = 0; lo < elements.size(); lo += wave_size)
    {
      const std::size_t hi = (std::min)(elements.size(), lo + wave_size);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(lo, hi),
        [&](const tbb::blocked_range<std::size_t>& r)
        {
          for (std::size_t i = r.begin(); i != r.end(); ++i)
          {
            const Element_type& element = elements[i];
            while (!op.lock_zone(element, c3t3))
            {
              c3t3.triangulation().unlock_all_elements();
              std::this_thread::yield();
            }
            op.execute_operation(element, c3t3);
            c3t3.triangulation().unlock_all_elements();
          }
        });
    }

    return true;
  }

  // Unordered processing (Vertex_smooth, Edge_flip): no ordering constraint.
  // Partition the candidate set into equal-count, spatially-compact buckets via
  // recursive median (kd-tree) splitting, then process buckets concurrently
  // with each bucket sequential on one thread. Equal counts -> balanced load
  // (no idle tail); compact boxes -> lock conflicts stay confined to bucket
  // boundaries (interiors never contend).
  bool apply_unordered_processing(std::vector<Element_type>& elements, Operation& op, C3t3& c3t3) const
  {
    CGAL_TR_SUB_VARS();
    CGAL_TR_SUB_T0();
    struct EP { double c[3]; Element_type e; };
    std::vector<EP> eps(elements.size());
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, elements.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t i = r.begin(); i != r.end(); ++i)
        {
          const auto p = op.point_on_element(elements[i]);
          eps[i] = EP{{CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())}, elements[i]};
        }
      });

    const std::size_t nthreads =
      static_cast<std::size_t>((std::max)(1, tbb::this_task_arena::max_concurrency()));
    const std::size_t target =
      (std::max)(std::size_t(1), eps.size() / (buckets_per_thread(op) * nthreads));

    constexpr std::size_t PAR_SPLIT_THRESHOLD = 8192;
    tbb::concurrent_vector<std::pair<std::size_t, std::size_t>> buckets;
    std::function<void(std::size_t, std::size_t)> split_range =
      [&](std::size_t lo, std::size_t hi)
      {
        if (hi - lo <= target)
        {
          if (hi > lo)
            buckets.emplace_back(lo, hi);
          return;
        }
        double mn[3] = {eps[lo].c[0], eps[lo].c[1], eps[lo].c[2]};
        double mx[3] = {eps[lo].c[0], eps[lo].c[1], eps[lo].c[2]};
        for (std::size_t i = lo + 1; i < hi; ++i)
          for (int a = 0; a < 3; ++a)
          {
            mn[a] = (std::min)(mn[a], eps[i].c[a]);
            mx[a] = (std::max)(mx[a], eps[i].c[a]);
          }
        int axis = 0;
        double best = mx[0] - mn[0];
        for (int a = 1; a < 3; ++a)
          if (mx[a] - mn[a] > best) { best = mx[a] - mn[a]; axis = a; }

        const std::size_t mid = lo + (hi - lo) / 2;
        std::nth_element(eps.begin() + lo, eps.begin() + mid, eps.begin() + hi,
                         [axis](const EP& x, const EP& y) { return x.c[axis] < y.c[axis]; });
        if (hi - lo > PAR_SPLIT_THRESHOLD)
        {
          tbb::task_group tg;
          tg.run([&, lo, mid] { split_range(lo, mid); });
          split_range(mid, hi);
          tg.wait();
        }
        else
        {
          split_range(lo, mid);
          split_range(mid, hi);
        }
      };
    split_range(0, eps.size());
    CGAL_TR_SUB_ADD(1); CGAL_TR_SUB_T0();

    // --- interior/boundary classification for lock elision (topological) ---
    // Tag each element-endpoint vertex with its bucket id (MIXED if >=2 buckets own
    // it); then flag every vertex incident to a cell whose tagged vertices span >=2
    // buckets ("boundary-touching"). An element is INTERIOR iff none of its locked
    // vertices is boundary-touching => its whole 1-ring lock zone is single-bucket, so
    // no other concurrently-processed bucket can touch it and it may run lock-free.
    // flip/smooth create/destroy no vertices, so this classification stays valid for
    // the whole phase. Empty locked_vertices() => the element falls through to the
    // locked path (the safe default).
    const std::size_t MIXED_BUCKET = static_cast<std::size_t>(-1);
    boost::concurrent_flat_map<Vertex_handle, std::size_t, boost::hash<Vertex_handle>> vertex_bucket;
    // Reserve to the true upper bound on distinct keys (every mesh vertex),
    // not a multiple of eps.size(): locked_vertices() now includes flip's ring
    // vertices, not just {v0, v1} (see Internal_edge_flip_operation in
    // flip_edges.h), so the same vertex is named by many different candidates
    // and an eps.size()-scaled guess undercounts badly. Sizing off eps.size()
    // alone left this map growing under heavy concurrent insertion, which
    // surfaced as std::bad_alloc from inside table_core -- confirmed by
    // backtrace on mesh3_243015 (see the "mesh3_243015 flip defect" log
    // entries) even after an *8 bump.
    vertex_bucket.reserve(c3t3.triangulation().number_of_vertices() + 64);

    // Pass 1: scatter element-endpoint vertex -> bucket id (MIXED on conflict).
    // locked_vertices() is cached per element here and reused by Pass 3 below
    // instead of being recomputed: for flip it circulates the edge's incident
    // facets, and computing it twice per candidate scanned (once here, once in
    // Pass 3) doubled that cost for no reason once the set stopped being the
    // trivial {v0, v1} pair.
    std::vector<boost::container::small_vector<Vertex_handle, 2>> all_lv(eps.size());
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, buckets.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t bi = r.begin(); bi != r.end(); ++bi)
          for (std::size_t i = buckets[bi].first; i < buckets[bi].second; ++i)
          {
            boost::container::small_vector<Vertex_handle, 2>& lv = all_lv[i];
            op.locked_vertices(eps[i].e, c3t3, lv);
            for (const Vertex_handle& v : lv)
              vertex_bucket.insert_or_visit(std::make_pair(v, bi),
                [bi, MIXED_BUCKET](std::pair<const Vertex_handle, std::size_t>& kv)
                { if (kv.second != bi) kv.second = MIXED_BUCKET; });
          }
      });

    CGAL_TR_SUB_ADD(2); CGAL_TR_SUB_T0();

    // Build the finite-cell list once (mesh is read-only during this preprocessing).
    std::vector<typename C3t3::Triangulation::Cell_handle> all_cells;
    {
      const auto& tr_ro = c3t3.triangulation();
      all_cells.reserve(tr_ro.number_of_finite_cells() + 64);
      for (auto cit = tr_ro.finite_cells_begin(); cit != tr_ro.finite_cells_end(); ++cit)
        all_cells.push_back(cit);
    }
    CGAL_TR_SUB_ADD(3); CGAL_TR_SUB_T0();

    // Pass 2: a cell whose tagged vertices span >=2 buckets (any MIXED, or two
    // different single-bucket tags) is shared; flag all 4 of its vertices
    // boundary-touching. Untagged vertices are ignored (own no element -> no lock).
    boost::concurrent_flat_map<Vertex_handle, char, boost::hash<Vertex_handle>> boundary_vertex;
    boundary_vertex.reserve(c3t3.triangulation().number_of_vertices() + 64);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, all_cells.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t ci = r.begin(); ci != r.end(); ++ci)
        {
          const auto c = all_cells[ci];
          std::size_t first = 0; bool have = false, shared = false;
          for (int k = 0; k < 4; ++k)
          {
            std::size_t b = 0;
            const bool tagged = vertex_bucket.cvisit(c->vertex(k),
              [&b](const std::pair<const Vertex_handle, std::size_t>& kv){ b = kv.second; }) > 0;
            if (!tagged) continue;
            if (b == MIXED_BUCKET || (have && b != first)) { shared = true; break; }
            first = b; have = true;
          }
          if (shared)
            for (int k = 0; k < 4; ++k)
              boundary_vertex.insert_or_assign(c->vertex(k), char(1));
        }
      });

    CGAL_TR_SUB_ADD(4); CGAL_TR_SUB_T0();

    // Pass 2b: grow the boundary tag by one cell-layer -- ONLY for operations whose
    // footprint reaches past their locked_vertices' 1-rings. Restored to a single
    // growth pass: flip's locked_vertices() now includes its edge's ring vertices
    // (see Internal_edge_flip_operation::locked_vertices() in flip_edges.h), so Pass 2
    // above already gives full-star coverage for any ring vertex actually shared
    // across buckets -- a cell containing a MIXED-tagged ring vertex is flagged
    // directly by Pass 2's own per-cell check, for every one of that vertex's
    // incident cells, without needing growth. What growth still needs to cover is
    // flip_n_to_m's neighbour-cell write, one hop beyond the edge-incident cells
    // (see the "mesh3_243015 flip defect" / "locked_vertices widened" log entries for
    // the full footprint audit) -- bounded, so one layer suffices. (An earlier
    // attempt widened this pass to a fixed point instead of widening
    // locked_vertices(); it was correct but cost ~2x the flip phase because the
    // fixed point of "shares a cell with a tagged vertex" swallows nearly the whole
    // mesh once real interfaces exist. Kept as a cautionary note, not code.)
    // Must happen here, before Pass 3, which concurrently mutates the TDS this reads.
    if (op.footprint_needs_halo())
    {
      boost::concurrent_flat_map<Vertex_handle, char, boost::hash<Vertex_handle>> halo;
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, all_cells.size()),
        [&](const tbb::blocked_range<std::size_t>& r)
        {
          for (std::size_t ci = r.begin(); ci != r.end(); ++ci)
          {
            const auto c = all_cells[ci];
            bool touches = false;
            for (int k = 0; k < 4 && !touches; ++k)
              if (boundary_vertex.contains(c->vertex(k)))
                touches = true;
            if (touches)
              for (int k = 0; k < 4; ++k)
                halo.insert_or_assign(c->vertex(k), char(1));
          }
        });
      halo.cvisit_all([&](const std::pair<const Vertex_handle, char>& kv)
                      { boundary_vertex.insert_or_assign(kv.first, char(1)); });
    }
    CGAL_TR_SUB_ADD(5); CGAL_TR_SUB_T0();

    // Pass 3: process buckets concurrently; each bucket sequential on its worker.
    // INTERIOR elements skip lock_zone (lock-free) and execute immediately.
    // BOUNDARY elements: if op.prefers_boundary_waves(), collect the index
    // instead of executing -- handled below, once every interior element in
    // this phase has finished, by a conflict-free-wave pass scoped to just
    // this (typically small) subset. Otherwise, unchanged: lock and execute
    // inline via op.lock_zone_precomputed()/lock_zone().
#ifdef CGAL_TR_SUBPHASE_TIMINGS
    std::atomic<std::size_t> n_int_{0}, n_lck_{0}, n_rty_{0};
#endif
    const bool boundary_waves = op.prefers_boundary_waves();
    const bool bucket_scoped_locks = op.prefers_bucket_scoped_locks();
    tbb::combinable<std::vector<std::size_t>> local_boundary_idx;
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, buckets.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t bi = r.begin(); bi != r.end(); ++bi)
        {
          for (std::size_t i = buckets[bi].first; i < buckets[bi].second; ++i)
          {
            const Element_type& element = eps[i].e;
            const boost::container::small_vector<Vertex_handle, 2>& lv = all_lv[i];
            bool interior = !lv.empty();
            for (const Vertex_handle& v : lv)
              if (boundary_vertex.contains(v)) { interior = false; break; }

            if (interior)
            {
              CGAL_TR_COUNT(interior);
#ifdef CGAL_TR_SUBPHASE_TIMINGS
              n_int_.fetch_add(1, std::memory_order_relaxed);
#endif
              CGAL_TR_PROBE_MODE(1);
              op.execute_operation(element, c3t3);
              // An operation may acquire vertex locks on demand even on the
              // lock-elided path (see flip_edges.h, tr_flip_require_lock):
              // release them here or they leak into the next element --
              // unless the operation asked to keep them for the whole bucket,
              // in which case they are released once the bucket is done.
#ifndef CGAL_TR_FLIP_NO_LOCK_ON_DEMAND
              if (!bucket_scoped_locks)
              {
                c3t3.triangulation().unlock_all_elements();
                op.on_locks_released();
              }
#endif
            }
            else if (boundary_waves)
            {
              CGAL_TR_COUNT(locked);
              local_boundary_idx.local().push_back(i);
            }
            else
            {
              CGAL_TR_COUNT(locked);
#ifdef CGAL_TR_SUBPHASE_TIMINGS
              n_lck_.fetch_add(1, std::memory_order_relaxed);
#endif
              while (!op.lock_zone_precomputed(element, c3t3, lv))
              {
                CGAL_TR_COUNT(lock_retry);
#ifdef CGAL_TR_SUBPHASE_TIMINGS
                n_rty_.fetch_add(1, std::memory_order_relaxed);
#endif
                c3t3.triangulation().unlock_all_elements();
                op.on_locks_released();
                std::this_thread::yield();
              }
              {
                CGAL_TR_PROBE_MODE(2);
                op.execute_operation(element, c3t3);
              }
              if (!bucket_scoped_locks)
              {
                c3t3.triangulation().unlock_all_elements();
                op.on_locks_released();
              }
            }
          }

          // End of bucket: release everything held across this bucket's
          // elements and drop whatever the operation cached under those locks.
          if (bucket_scoped_locks)
          {
            c3t3.triangulation().unlock_all_elements();
            op.on_locks_released();
          }
        }
      });

    CGAL_TR_SUB_ADD(6); CGAL_TR_SUB_T0();

    // Boundary conflict-free-wave pass (only entered when
    // op.prefers_boundary_waves() is true): same wave-selection algorithm as
    // apply_conflict_free_waves() above, but over just the deferred boundary
    // subset, and reusing the Pass-1-cached all_lv instead of recomputing
    // locked_vertices() per element. See prefers_boundary_waves()'s comment in
    // the base class for the soundness argument (the interior/boundary split
    // already guarantees this subset's footprint cannot have been touched by
    // the interior work that just finished above).
    if (boundary_waves)
    {
      std::vector<std::size_t> pending;
      local_boundary_idx.combine_each([&](const std::vector<std::size_t>& v)
        { pending.insert(pending.end(), v.begin(), v.end()); });

      if (!pending.empty())
      {
        using Cell_handle = typename C3t3::Triangulation::Cell_handle;
        auto& tr = c3t3.triangulation();
        std::vector<std::size_t> wave, deferred;
        boost::container::small_vector<Cell_handle, 64> inc_cells;
        boost::container::small_vector<Vertex_handle, 16> ext;
        const bool needs_halo = op.footprint_needs_halo();

        while (!pending.empty())
        {
          const std::size_t wave_stamp = next_wave_claim_stamp();
          wave.clear();
          deferred.clear();

          for (const std::size_t idx : pending)
          {
            const boost::container::small_vector<Vertex_handle, 2>& lv = all_lv[idx];
            if (lv.empty()) { deferred.push_back(idx); continue; }

            // The wave claim/test must cover the SAME footprint Pass 2b's
            // halo growth covers for classification (footprint_needs_halo()):
            // flip's actual read/write reaches one cell layer beyond
            // locked_vertices()'s own star, not just the star itself. Testing
            // and stamping only lv's 1-ring (as apply_conflict_free_waves()
            // does for operations that don't need a halo) missed exactly this
            // extra hop -- confirmed by a SIGSEGV/livelock on the first
            // attempt at this scheme that stamped only lv's 1-ring.
            //
            // Growing from EVERY vertex in lv (including lv[0]/lv[1], the
            // element's own two endpoints for the one operation that uses
            // this path today) was tried and measured as a severe
            // over-claim: an edge endpoint's incident-cell star reaches
            // every cell touching that vertex anywhere in the mesh, not just
            // the cells near this specific edge, so two elements sharing a
            // busy, unrelated endpoint would spuriously "conflict" even with
            // disjoint true footprints -- confirmed by 13/13 runs timing out
            // (see the "mesh3_243015 flip defect" / ATTEMPT 9 log entry).
            // Pass 2b's classification halo does not have this problem
            // because it grows from CELLS already known to be near an
            // interface, not from a vertex's whole star.
            //
            // Contract (documented here since it is a real coupling, not
            // enforced by the type system): an operation using
            // prefers_boundary_waves() together with footprint_needs_halo()
            // must order locked_vertices()'s output as [own primary
            // vertices][ring/secondary vertices] and pass how many leading
            // entries are primary via halo_growth_skip() (default 0, meaning
            // "grow from all of lv" -- correct for a widened claim whose
            // vertices are all inherently local, just not for endpoints like
            // flip's e.first/e.second whose stars are unbounded).
            ext.clear();
            for (const Vertex_handle& v : lv)
              if (std::find(ext.begin(), ext.end(), v) == ext.end())
                ext.push_back(v);
            if (needs_halo)
            {
              const std::size_t skip = (std::min)(op.halo_growth_skip(), lv.size());
              for (std::size_t li = skip; li < lv.size(); ++li)
              {
                const Vertex_handle& v = lv[li];
                inc_cells.clear();
                tr.incident_cells(v, std::back_inserter(inc_cells));
                for (const Cell_handle& c : inc_cells)
                  for (int k = 0; k < 4; ++k)
                  {
                    const Vertex_handle w = c->vertex(k);
                    if (std::find(ext.begin(), ext.end(), w) == ext.end())
                      ext.push_back(w);
                  }
              }
            }

            bool conflicts = false;
            for (const Vertex_handle& v : ext)
              if (v->wave_claim_stamp() == wave_stamp) { conflicts = true; break; }
            if (conflicts) { deferred.push_back(idx); continue; }

            for (const Vertex_handle& v : ext)
            {
              inc_cells.clear();
              tr.incident_cells(v, std::back_inserter(inc_cells));
              for (const Cell_handle& c : inc_cells)
                for (int k = 0; k < 4; ++k)
                  c->vertex(k)->set_wave_claim_stamp(wave_stamp);
            }
            wave.push_back(idx);
          }

          if (wave.empty())
          {
            for (const std::size_t idx : deferred)
              op.execute_operation(eps[idx].e, c3t3);
            break;
          }

          // Sequential, not tbb::parallel_for: a crashbt-instrumented run
          // caught a genuine crash inside find_best_flip()'s lazy population
          // of the per-worker inc_cells_map cache while wave members ran in
          // parallel here -- meaning find_best_flip's real read footprint
          // (which explores/caches vertices while evaluating candidate flip
          // rotations) reaches further than the halo-extended `ext` claim
          // above accounts for, so two wave members that the claim considered
          // conflict-free could still, in practice, read/write overlapping
          // triangulation state. This only serializes the WAVE's own accepted
          // members (already a scoped, typically small minority of all
          // candidates in the phase -- the boundary subset only) rather than
          // falling back to sequential execution for the whole operation or
          // phase; interior elements above still run fully in parallel.
          //
          // A follow-up diagnostic (ATTEMPT 10, see flip_edges.h's
          // Internal_edge_flip_operation) tried tbb::parallel_for here again
          // with a scoped halo_growth_skip() fix applied first, to check
          // whether the earlier crash was really about over-broad growth
          // rather than parallel execution per se. It crashed identically --
          // confirming the gap is NOT about growth being too narrow or too
          // broad (scoped or unscoped both crash under parallel execution;
          // both are clean under sequential), i.e. no vertex-based claim this
          // scheduling layer can compute closes it. Sequential is correctly
          // the only validated-safe choice for this path.
          for (const std::size_t idx : wave)
            op.execute_operation(eps[idx].e, c3t3);

          pending.swap(deferred);
        }
      }
    }

    CGAL_TR_SUB_ADD(7);
#ifdef CGAL_TR_SUBPHASE_TIMINGS
    std::cout << "[sub] " << op.operation_name()
              << " | partition " << sub_[1]
              << " | pass1 "     << sub_[2]
              << " | celllist "  << sub_[3]
              << " | pass2 "     << sub_[4]
              << " | pass2b "    << sub_[5]
              << " | pass3 "     << sub_[6]
              << " | bwaves "    << sub_[7]
              << " | buckets "   << buckets.size()
              << " | elements "  << eps.size()
              << " | interior "  << n_int_.load()
              << " | locked "    << n_lck_.load()
              << " | retries "   << n_rty_.load()
              << std::endl;
#endif
    return true;
  }
};

#endif // CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && CGAL_LINKED_WITH_TBB

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H
