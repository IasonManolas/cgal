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

#include <CGAL/tags.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/number_utils.h>
#include <unordered_set>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/Tetrahedral_remeshing/internal/Parallel_tuning.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/blocked_range.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_unordered_map.h>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered/concurrent_flat_map.hpp>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/task_arena.h>
#endif

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdio>
#include <chrono>
#include <numeric>
#include <dlfcn.h>
#include <iterator>
#include <string>
#include <thread>
#include <cmath>
#include <iterator>
#include <cstdlib>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <vector>

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#include <CGAL/Real_timer.h>
#include <cstddef>
#include <iostream>
#endif

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

#ifdef CGAL_LINKED_WITH_TBB
/**
* D2 -- one spatial subdivision per remeshing iteration, shared by every
* unordered operation in it.
*
* run_unordered() currently re-partitions per operation: internal flip,
* boundary flip and three smoothing passes each run their own nth_element
* recursion over their own candidate set, five times an iteration. The element
* TYPES differ (vertex pairs, single vertices) but the geometry does not --
* all of them are points in the same triangulation, and what is being computed
* is a spatial subdivision of it.
*
* So compute the subdivision once and keep only its split PLANES. Later
* operations assign their candidates by walking the tree, Theta(candidates)
* with no recursion and no median selection. Two effects, and the second is
* probably the larger: the partition cost is paid once instead of five times,
* and every operation in the iteration uses THE SAME regions, so a thread that
* worked on a region during flip works on it again during smoothing with that
* part of the mesh still in its cache. There is no cross-phase locality today.
*
* A perfect binary tree, so the leaf count is a power of two.
*/
struct Shared_kd_partition
{
  struct Node { int axis = 0; double value = 0.; };
  std::vector<Node> nodes;      // implicit tree: children of i are 2i+1, 2i+2
  std::size_t leaves = 0;
  bool valid = false;

  void reset() { nodes.clear(); leaves = 0; valid = false; }

  std::size_t leaf_of(const double x, const double y, const double z) const
  {
    std::size_t i = 0;
    const std::size_t internal = leaves - 1;
    while (i < internal)
    {
      const Node& n = nodes[i];
      const double c = (n.axis == 0) ? x : (n.axis == 1) ? y : z;
      i = (c < n.value) ? (2 * i + 1) : (2 * i + 2);
    }
    return i - internal;
  }

  static Shared_kd_partition& get()
  {
    static Shared_kd_partition s;
    return s;
  }
};

/**
* D1 -- METIS, loaded at run time.
*
* dlopen rather than a link-time dependency on purpose: the measurement build
* must not change shape between the arms (POLICY 0.2), and CGAL's
* Tetrahedral_remeshing gaining a hard external dependency is a packaging
* decision that should follow a measured result, not precede it. If the library
* is absent or the call fails, the candidate falls back to the shipped kd
* partition and says so through the return value, so a missing METIS shows up
* as "no change", never as a wrong partition.
*/
struct Metis_api
{
  using PartKway = int (*)(int*, int*, int*, int*, int*, int*, int*, int*,
                           float*, float*, int*, int*, int*);
  using SetOpts  = int (*)(int*);
  PartKway part_kway = nullptr;
  SetOpts  set_opts  = nullptr;
  bool ok = false;

  static const Metis_api& get()
  {
    static const Metis_api api = load();
    return api;
  }

private:
  static Metis_api load()
  {
    Metis_api a;
    void* h = dlopen("libmetis.so.5", RTLD_LAZY);
    if (h == nullptr)
      h = dlopen("libmetis.so", RTLD_LAZY);
    if (h == nullptr)
      return a;
    a.part_kway = reinterpret_cast<PartKway>(dlsym(h, "METIS_PartGraphKway"));
    a.set_opts  = reinterpret_cast<SetOpts>(dlsym(h, "METIS_SetDefaultOptions"));
    a.ok = (a.part_kway != nullptr && a.set_opts != nullptr);
    return a;
  }
};
#endif // CGAL_LINKED_WITH_TBB

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
* elements first. `false` groups them by a coarse spatial grid instead, so
* that threads work in different regions of the triangulation and lock
* conflicts stay down.
*/
#ifdef CGAL_LINKED_WITH_TBB
template <typename Operation>
class Elementary_operation_execution_parallel
{
public:
  using C3t3 = typename Operation::C3t3;
  using Element_type = typename Operation::Element_type;
  using Element_range = typename Operation::Element_range;

  bool execute(Operation& op, C3t3& c3t3) const
  {
    std::vector<Element_type> candidates = collect(op, c3t3);
    if (candidates.empty())
      return false;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    CGAL::Real_timer timer;
    timer.start();
    const std::size_t nb_candidates = candidates.size();
#endif

    if constexpr (Operation::requires_ordered_processing)
      run_ordered(candidates, op, c3t3);
    else
      run_unordered(candidates, op, c3t3);

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
#ifdef CGAL_TR_ZONE_STATS
    ++Zone_stats::get().zones_locked;
    ++Zone_stats::get().zone_attempts;
    // Time only every 64th failed attempt. Reading a clock on a path taken
    // millions of times per run would perturb what it measures; counts stay
    // exact, the nanoseconds are sampled.
    std::size_t zs_fail = 0;
    const bool zs_time = ((Zone_stats::get().zone_attempts.load() & 63u) == 0u);
    const auto zs_t0 = std::chrono::steady_clock::now();
#endif
    while (!op.lock_zone(element, c3t3))
    {
#ifdef CGAL_TR_ZONE_STATS
      ++Zone_stats::get().zone_attempts;
      ++Zone_stats::get().yields;
      ++zs_fail;
#endif
      c3t3.triangulation().unlock_all_elements();
      std::this_thread::yield();
    }
#ifdef CGAL_TR_ZONE_STATS
    if (zs_time && zs_fail)
    {
      Zone_stats::get().ns_sampled += static_cast<std::size_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - zs_t0).count());
      Zone_stats::get().n_sampled += zs_fail;
    }
#endif
    op.execute_operation(element, c3t3);
    c3t3.triangulation().unlock_all_elements();
  }

  /**
  * `CGAL_TR_DEFER_ON_CONFLICT=1` (A6). Two problems with the spin above, and
  * the second is the expensive one:
  *
  *  1. `yield()` returns immediately and retries, so a thread waiting on a
  *     contended zone burns a core. That inflates task-clock without doing
  *     work, which means the measured CPU utilisation is flattering and the
  *     real idle fraction is larger than it looks.
  *  2. Every retry THROWS AWAY the star walk. `lock_zone()` calls
  *     `try_lock_and_get_incident_cells()` twice, which BFS-walks both stars
  *     and marks and unmarks `tds_data()`, and then discards all of it when a
  *     later lock in the same zone fails. Under contention the same star is
  *     walked over and over.
  *
  * So: try a bounded number of times, then put the element aside and take the
  * next one. The bucket runs its deferred list at the end, with the unbounded
  * loop as a backstop, so no element is ever dropped and termination is
  * unchanged.
  *
  * Only the UNORDERED path may do this. Deferring changes the order elements
  * are processed in, and split depends on the global longest-first order --
  * R2 lost 6.4% when that was broken. Flip and smooth have no such dependency;
  * `run_ordered()` and the collapse executor are untouched.
  */
  static void apply_one_or_defer(const Element_type& element, Operation& op,
                                 C3t3& c3t3,
                                 std::vector<Element_type>& deferred)
  {
    for (int attempt = 0; attempt < 2; ++attempt)
    {
      if (op.lock_zone(element, c3t3))
      {
        op.execute_operation(element, c3t3);
        c3t3.triangulation().unlock_all_elements();
        return;
      }
      c3t3.triangulation().unlock_all_elements();
    }
    deferred.push_back(element);
  }

  /**
  * `CGAL_TR_BUCKET_ORDERED=1` groups the ordered elements too. It is off by
  * default and is a different trade from the unordered case: draining one
  * queue keeps the order `get_elements()` produced -- longest edge first for
  * split -- across all threads, whereas grouping keeps it only WITHIN a
  * bucket, and buckets run concurrently. So this buys the same locality at
  * the cost of the global ordering, and only measurement says whether that is
  * worth it.
  *
  * Only split reaches this. Collapse has its own executor: its work list
  * changes as it runs and is drained from a priority queue that the collapses
  * themselves push back into, which grouping cannot express.
  */
  static bool kd_ordered_enabled()
  {
    static const bool enabled = []
      {
        const char* const e = std::getenv("CGAL_TR_KD_ORDERED");
        return (e != nullptr) && (std::atoi(e) != 0);
      }();
    return enabled;
  }

  static bool bucket_ordered_enabled()
  {
    static const bool enabled = []
      {
        const char* const e = std::getenv("CGAL_TR_BUCKET_ORDERED");
        return (e != nullptr) && (std::atoi(e) != 0);
      }();
    return enabled;
  }

  static void run_ordered(std::vector<Element_type>& candidates,
                          Operation& op, C3t3& c3t3)
  {
    // R2 measured grouping the ordered elements by a uniform grid: -7.170%,
    // because grouping keeps the order only WITHIN a bucket and split depends
    // on the global longest-first order. `CGAL_TR_KD_ORDERED=1` asks whether
    // equal-count buckets change that answer -- they fixed the unordered case
    // -- or whether the loss is the lost ordering rather than the lopsided
    // buckets. Both off by default.
    if (kd_ordered_enabled())
    {
      std::vector<std::vector<Element_type> > parts = kd_partition(candidates, op);
      return run_parts(parts, op, c3t3);
    }

    if (bucket_ordered_enabled())
    {
      Buckets buckets = bucket_by_grid(candidates, op, c3t3);
      if (!buckets.empty())
        return run_buckets(buckets, op, c3t3);
    }

    tbb::concurrent_queue<Element_type> queue(candidates.begin(), candidates.end());
    tbb::parallel_for(0, tbb::this_task_arena::max_concurrency(),
                      [&](int)
                      {
                        Element_type element;
                        while (queue.try_pop(element))
                          apply_one(element, op, c3t3);
                      });
  }

  struct Grid_cell_index
  {
    int i, j, k;
    bool operator==(const Grid_cell_index& o) const
    { return i == o.i && j == o.j && k == o.k; }
  };
  struct Grid_cell_index_hash
  {
    std::size_t operator()(const Grid_cell_index& ci) const
    {
      return ((std::hash<int>()(ci.i) ^ (std::hash<int>()(ci.j) << 1)) >> 1)
           ^ (std::hash<int>()(ci.k) << 1);
    }
  };
  using Buckets = std::unordered_map<Grid_cell_index, std::vector<Element_type>,
                                     Grid_cell_index_hash>;

  /**
  * Groups the elements by the cell of a coarse grid, keeping each bucket in
  * the order the elements arrived. Returns an empty map if the bounding box
  * is degenerate, which the callers read as "do not group".
  */
  static Buckets bucket_by_grid(const std::vector<Element_type>& candidates,
                                const Operation& op, const C3t3& c3t3)
  {
    const CGAL::Bbox_3 bb = c3t3.bbox();
    const double min_sq_dim
      = (std::min)(CGAL::square(bb.xmax() - bb.xmin()),
        (std::min)(CGAL::square(bb.ymax() - bb.ymin()),
                   CGAL::square(bb.zmax() - bb.zmin())));
    const double cell_size = 0.5 * CGAL::approximate_sqrt(min_sq_dim);

    Buckets buckets;
    if (cell_size <= 0.)
      return buckets;

    const double inv_cell_size = 1. / cell_size;
    for (const Element_type& element : candidates)
    {
      const auto p = op.point_on_element(element);
      buckets[Grid_cell_index{
                static_cast<int>(std::floor(p.x() * inv_cell_size)),
                static_cast<int>(std::floor(p.y() * inv_cell_size)),
                static_cast<int>(std::floor(p.z() * inv_cell_size))}]
        .push_back(element);
    }
    return buckets;
  }

  /**
  * Splits `range` at the median of its widest extent, recursively, until there
  * are at least `target` parts. Each part is appended to `out`.
  */
  static void kd_split(typename std::vector<Element_type>::iterator first,
                       typename std::vector<Element_type>::iterator last,
                       const Operation& op, std::size_t target,
                       std::vector<std::vector<Element_type> >& out)
  {
    const std::size_t n = static_cast<std::size_t>(std::distance(first, last));
    if (n == 0)
      return;
    if (target <= 1 || n < 2)
    {
      out.emplace_back(first, last);
      return;
    }

    // widest coordinate extent decides the split axis
    double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
    double hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
    for (auto it = first; it != last; ++it)
    {
      const auto p = op.point_on_element(*it);
      const double c[3] = { CGAL::to_double(p.x()),
                            CGAL::to_double(p.y()),
                            CGAL::to_double(p.z()) };
      for (int a = 0; a < 3; ++a)
      {
        lo[a] = (std::min)(lo[a], c[a]);
        hi[a] = (std::max)(hi[a], c[a]);
      }
    }
    int axis = 0;
    for (int a = 1; a < 3; ++a)
      if (hi[a] - lo[a] > hi[axis] - lo[axis])
        axis = a;

    if (!(hi[axis] > lo[axis])) // all in one place: nothing to split by
    {
      out.emplace_back(first, last);
      return;
    }

    const auto coord = [&](const Element_type& e)
    {
      const auto p = op.point_on_element(e);
      return (axis == 0) ? CGAL::to_double(p.x())
           : (axis == 1) ? CGAL::to_double(p.y())
                         : CGAL::to_double(p.z());
    };

    const auto mid = first + static_cast<std::ptrdiff_t>(n / 2);
    std::nth_element(first, mid, last,
                     [&](const Element_type& a, const Element_type& b)
                     { return coord(a) < coord(b); });

    kd_split(first, mid, op, target / 2, out);
    kd_split(mid, last, op, target - target / 2, out);
  }

  static std::vector<std::vector<Element_type> >
  kd_partition(std::vector<Element_type>& candidates, const Operation& op)
  {
    // four buckets per thread, so work-stealing has something to steal.
    // D4a: `CGAL_TR_BUCKETS_PER_THREAD` makes that a measured choice rather
    // than a guess. It sets the TAIL: a phase ends when its last bucket ends,
    // so with b buckets per thread the tail is about 1/b of a thread's share.
    // More buckets shorten it, and cost more partitioning, more interface
    // between buckets, and less elidable interior.
    const std::size_t target =
      static_cast<std::size_t>(Parallel_tuning::get().buckets_per_thread)
      * static_cast<std::size_t>((std::max)(1, tbb::this_task_arena::max_concurrency()));

    std::vector<std::vector<Element_type> > parts;
    parts.reserve(target);
    kd_split(candidates.begin(), candidates.end(), op, target, parts);
    return parts;
  }

  /**
  * `CGAL_TR_ELIDE_INTERIOR_LOCKS=1` lets an element whose lock zone cannot be
  * reached by any other bucket run without taking any lock. Off by default.
  *
  * The classification is topological, not geometric. Tag each element's locked
  * vertices with the bucket that owns them, MIXED if two buckets do. A cell
  * whose tagged vertices span more than one bucket is shared, so flag all four
  * of its vertices as boundary-touching. An element is interior when none of
  * its locked vertices is flagged: every cell around those vertices then
  * belongs to one bucket, and since a bucket is one task, no other thread can
  * be inside that zone.
  *
  * The argument covers the ONE-RING of the locked vertices, and it is only
  * valid for an operation whose write footprint is that one-ring. Smoothing
  * qualifies: it writes the vertex and re-checks the orientation of the cells
  * incident to it. Flip does NOT -- it re-stitches the mirror cells across the
  * star's outer facets, which live in the two-ring, which is exactly why
  * lock_flip_zone() locks the star's neighbours as well. An operation that
  * does not opt in returns no locked vertices and always takes its locks.
  */
  static bool elide_interior_locks()
  {
    static const bool enabled = []
      {
        const char* const e = std::getenv("CGAL_TR_ELIDE_INTERIOR_LOCKS");
        return (e != nullptr) && (std::atoi(e) != 0);
      }();
    return enabled;
  }


  /**
  * CGAL_TR_ELISION_STATS=1 reports, per operation and per call, how many
  * elements the interior/boundary classification actually freed from their
  * lock zone. It exists because the classification is paid unconditionally
  * (the elision_mode switch falls through to it), and its cost is only
  * justified by what it produces. halofix measured `interior 0 | locked
  * 392764` on every internal-flip call -- the chain running dead. Whether the
  * same holds here is a measurement, not an inference, and this is it.
  */
  static bool elision_stats()
  {
    static const bool on = []{
      const char* const e = std::getenv("CGAL_TR_ELISION_STATS");
      return e != nullptr && *e == '1';
    }();
    return on;
  }

  using Vertex_handle = typename C3t3::Triangulation::Vertex_handle;
  using Interior_set = boost::concurrent_flat_map<Vertex_handle, char,
                                                  boost::hash<Vertex_handle> >;

  /**
  * Returns the vertices that are safe to work around without locking, or an
  * empty map if the operation does not support elision or it is switched off.
  */
  static Interior_set classify_interior(
      const std::vector<std::vector<Element_type> >& parts,
      const Operation& op, const C3t3& c3t3)
  {
    Interior_set interior;
    if (!elide_interior_locks())
      return interior;

    boost::container::small_vector<Vertex_handle, 2> probe;
    op.locked_vertices(parts.empty() || parts[0].empty()
                       ? Element_type() : parts[0][0], probe);
    if (probe.empty())
      return interior; // operation does not opt in

    static const std::size_t MIXED = static_cast<std::size_t>(-1);
    boost::concurrent_flat_map<Vertex_handle, std::size_t,
                               boost::hash<Vertex_handle> > owner;

    // which bucket owns each locked vertex, MIXED if more than one does
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, parts.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        boost::container::small_vector<Vertex_handle, 2> lv;
        for (std::size_t b = r.begin(); b != r.end(); ++b)
          for (const Element_type& e : parts[b])
          {
            lv.clear();
            op.locked_vertices(e, lv);
            for (const Vertex_handle& v : lv)
              owner.insert_or_visit(std::make_pair(v, b),
                [b](std::pair<const Vertex_handle, std::size_t>& kv)
                { if (kv.second != b) kv.second = MIXED; });
          }
      });

    // a cell spanning more than one bucket makes all four of its vertices
    // boundary-touching; untagged vertices own no element and lock nothing
    std::vector<typename C3t3::Triangulation::Cell_handle> cells;
    const auto& tr = c3t3.triangulation();
    cells.reserve(tr.number_of_finite_cells() + 64);
    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
      cells.push_back(cit);

    boost::concurrent_flat_map<Vertex_handle, char,
                               boost::hash<Vertex_handle> > touching;
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, cells.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t ci = r.begin(); ci != r.end(); ++ci)
        {
          const auto c = cells[ci];
          std::size_t first = 0;
          bool have = false, shared = false;
          for (int k = 0; k < 4; ++k)
          {
            std::size_t b = 0;
            const bool tagged = owner.cvisit(c->vertex(k),
              [&b](const std::pair<const Vertex_handle, std::size_t>& kv)
              { b = kv.second; }) > 0;
            if (!tagged)
              continue;
            if (b == MIXED || (have && b != first)) { shared = true; break; }
            first = b; have = true;
          }
          if (shared)
            for (int k = 0; k < 4; ++k)
              touching.insert_or_assign(c->vertex(k), char(1));
        }
      });

    owner.cvisit_all([&](const std::pair<const Vertex_handle, std::size_t>& kv)
    {
      if (kv.second != MIXED && !touching.contains(kv.first))
        interior.insert_or_assign(kv.first, char(1));
    });
    return interior;
  }

  /** True when every vertex this element locks is bucket-interior. */
  static bool is_interior(const Element_type& e, const Operation& op,
                          const Interior_set& interior)
  {
    if (interior.empty())
      return false;
    boost::container::small_vector<Vertex_handle, 2> lv;
    op.locked_vertices(e, lv);
    if (lv.empty())
      return false;
    for (const Vertex_handle& v : lv)
      if (!interior.contains(v))
        return false;
    return true;
  }

  /**
  * Runs the element without taking its lock zone. Only for an element the
  * classification says no other bucket can reach.
  */
  static void apply_one_unlocked(const Element_type& element,
                                 Operation& op, C3t3& c3t3)
  {
    op.execute_operation(element, c3t3);

    // The ZONE lock is what is elided, not every lock. The operation still
    // takes locks of its own -- smoothing locks the point it moves the vertex
    // to, and each intermediate position it tries on the way -- and those are
    // released here exactly as on the locked path. Without this a thread keeps
    // every grid cell it ever touched, and once it holds enough of them the
    // others cannot make progress: measured as a 6.2 second run still going
    // after 17 minutes.
    c3t3.triangulation().unlock_all_elements();
  }

  /**
  * The zone an element actually excludes other threads from, expressed as a
  * set of LOCK-GRID CELL indices (candidates B2..B4).
  *
  * Why grid cells and not vertices. R19 reasoned about vertices and one-rings,
  * and then had to argue separately for each operation about whether its write
  * footprint fits inside a one-ring -- which is why flip had to opt out
  * entirely, its footprint being the two-ring. The grid cell is the unit of
  * mutual exclusion itself: if no other bucket can touch any grid cell this
  * element's zone locks, no other thread can be inside the zone, whatever the
  * footprint is. `zone_ring` says how far the operation's footprint reaches
  * (1 for smoothing, 2 for flip), so the same code covers both.
  *
  * Conservative in the safe direction only: a coarse grid can call a genuinely
  * interior element shared, never the reverse.
  */
  template <typename OutVec>
  static void zone_grid_indices(const Element_type& e, const Operation& op,
                                const C3t3& c3t3, OutVec& out)
  {
    out.clear();
    boost::container::small_vector<Vertex_handle, 2> lv;
    op.elision_vertices(e, lv);
    if (lv.empty())
      return;

    using Cell_handle = typename C3t3::Triangulation::Cell_handle;
    const auto& tr = c3t3.triangulation();
    boost::container::small_vector<Cell_handle, 64> ring;
    for (const Vertex_handle& v : lv)
      tr.incident_cells(v, std::back_inserter(ring));

    const std::size_t one_ring = ring.size();
    if (Operation::zone_ring >= 2)
      for (std::size_t i = 0; i < one_ring; ++i)
        for (int k = 0; k < 4; ++k)
          ring.push_back(ring[i]->neighbor(k));

    for (const Cell_handle c : ring)
      for (int k = 0; k < 4; ++k)
      {
        const int gi = tr.lock_grid_index(c->vertex(k)->point());
        if (gi >= 0 && std::find(out.begin(), out.end(), gi) == out.end())
          out.push_back(gi);
      }
  }

  using Interior_flags = std::vector<std::vector<char> >;

  /**
  * B1 -- the same property R19 elides on, established WITHOUT a full-mesh scan.
  *
  * R19's classify_interior() walked every finite cell of the triangulation
  * into a std::vector and then ran two concurrent hash maps over it, once per
  * operation per phase. That is Theta(mesh) work, twice, to save locks on
  * Theta(candidates) elements -- and late in a remeshing run the candidates
  * are few and the mesh is not, which is exactly when the trade is worst. It
  * measured -0.98%.
  *
  * The predicate it was computing can be tested directly on each element: an
  * element of bucket b is interior when, for every vertex it locks, every
  * vertex of every cell in the relevant ring is either untagged (owns no
  * element, so locks nothing) or owned by b. That is the property the
  * independent elision validator checks, tested here rather than derived
  * through a `touching` map. Cost is Theta(sum of candidate rings); no cell
  * outside a candidate's ring is ever visited.
  */
  static Interior_flags classify_flags_vertex(
      const std::vector<std::vector<Element_type> >& parts,
      const Operation& op, const C3t3& c3t3)
  {
    static const std::size_t MIXED = static_cast<std::size_t>(-1);
    boost::concurrent_flat_map<Vertex_handle, std::size_t,
                               boost::hash<Vertex_handle> > owner;

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, parts.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        boost::container::small_vector<Vertex_handle, 2> lv;
        for (std::size_t b = r.begin(); b != r.end(); ++b)
          for (const Element_type& e : parts[b])
          {
            lv.clear();
            op.elision_vertices(e, lv);
            for (const Vertex_handle& v : lv)
              owner.insert_or_visit(std::make_pair(v, b),
                [b](std::pair<const Vertex_handle, std::size_t>& kv)
                { if (kv.second != b) kv.second = MIXED; });
          }
      });

    Interior_flags flags(parts.size());
    using Cell_handle = typename C3t3::Triangulation::Cell_handle;
    const auto& tr = c3t3.triangulation();

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, parts.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        boost::container::small_vector<Vertex_handle, 2> lv;
        boost::container::small_vector<Cell_handle, 64> ring;
        for (std::size_t b = r.begin(); b != r.end(); ++b)
        {
          flags[b].assign(parts[b].size(), char(0));
          for (std::size_t i = 0; i < parts[b].size(); ++i)
          {
            lv.clear();
            op.elision_vertices(parts[b][i], lv);
            if (lv.empty())
              continue;
            ring.clear();
            for (const Vertex_handle& v : lv)
              tr.incident_cells(v, std::back_inserter(ring));
            const std::size_t one_ring = ring.size();
            if (Operation::zone_ring >= 2)
              for (std::size_t j = 0; j < one_ring; ++j)
                for (int k = 0; k < 4; ++k)
                  ring.push_back(ring[j]->neighbor(k));

            bool interior = true;
            for (const Cell_handle c : ring)
            {
              for (int k = 0; k < 4 && interior; ++k)
              {
                std::size_t ob = 0;
                const bool tagged = owner.cvisit(c->vertex(k),
                  [&ob](const std::pair<const Vertex_handle, std::size_t>& kv)
                  { ob = kv.second; }) > 0;
                if (tagged && ob != b)
                  interior = false;
              }
              if (!interior)
                break;
            }
            flags[b][i] = interior ? char(1) : char(0);
          }
        }
      });
    return flags;
  }

  /**
  * B2 -- ownership of LOCK-GRID CELLS rather than of vertices.
  *
  * The zone index sets are computed ONCE and kept: the ownership pass and the
  * classification pass both need them, and recomputing them means walking
  * every candidate's ring twice. That doubling was measured (10.7 s against a
  * 4.5 s baseline on 118287_0.5) before the sets were cached.
  */
  static Interior_flags classify_flags_grid(
      const std::vector<std::vector<Element_type> >& parts,
      const Operation& op, const C3t3& c3t3)
  {
    static const std::size_t MIXED = static_cast<std::size_t>(-1);
    boost::concurrent_flat_map<int, std::size_t> owner;

    std::vector<std::vector<boost::container::small_vector<int, 16> > >
      zones(parts.size());

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, parts.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t b = r.begin(); b != r.end(); ++b)
        {
          zones[b].resize(parts[b].size());
          for (std::size_t i = 0; i < parts[b].size(); ++i)
          {
            zone_grid_indices(parts[b][i], op, c3t3, zones[b][i]);
            for (const int gi : zones[b][i])
              owner.insert_or_visit(std::make_pair(gi, b),
                [b](std::pair<const int, std::size_t>& kv)
                { if (kv.second != b) kv.second = MIXED; });
          }
        }
      });

    Interior_flags flags(parts.size());
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, parts.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t b = r.begin(); b != r.end(); ++b)
        {
          flags[b].assign(parts[b].size(), char(0));
          for (std::size_t i = 0; i < parts[b].size(); ++i)
          {
            const auto& idx = zones[b][i];
            if (idx.empty())
              continue;
            bool interior = true;
            for (const int gi : idx)
            {
              std::size_t ob = MIXED;
              owner.cvisit(gi, [&ob](const std::pair<const int, std::size_t>& kv)
                               { ob = kv.second; });
              if (ob != b) { interior = false; break; }
            }
            flags[b][i] = interior ? char(1) : char(0);
          }
        }
      });
    return flags;
  }

  static void run_flagged(std::vector<std::vector<Element_type> >& parts,
                          const Interior_flags& flags,
                          Operation& op, C3t3& c3t3)
  {
    const bool defer = Parallel_tuning::get().defer_on_conflict;
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, parts.size(), 1),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        std::vector<Element_type> deferred;
        for (std::size_t b = r.begin(); b != r.end(); ++b)
        {
          deferred.clear();
          for (std::size_t i = 0; i < parts[b].size(); ++i)
          {
            const Element_type& e = parts[b][i];
            if (!flags[b].empty() && flags[b][i])
              apply_one_unlocked(e, op, c3t3);
            else if (defer)
              apply_one_or_defer(e, op, c3t3, deferred);
            else
              apply_one(e, op, c3t3);
          }
          for (const Element_type& e : deferred)
            apply_one(e, op, c3t3);
        }
      });
  }

  /**
  * B3 -- acquire a bucket's whole grid-cell set once, then run every element
  * in it with no locking at all.
  *
  * A bucket holds candidates/(4*threads) elements, often hundreds, so paying
  * one acquisition instead of hundreds is the largest reduction in lock
  * traffic short of not locking. Ascending index order makes the multi-lock
  * acquisition deadlock-free without any global protocol.
  *
  * The failure mode is real and is R17's: if a bucket's grid-cell union is
  * most of the grid, no two buckets can run concurrently and the phase
  * serializes. Two guards, both required -- a cap on the union's size, and the
  * ordinary per-element path as a fallback that costs one failed attempt.
  */
  static void run_bucket_scoped(std::vector<std::vector<Element_type> >& parts,
                                Operation& op, C3t3& c3t3)
  {
    const std::size_t nb = parts.size();
    std::vector<std::vector<int> > unions(nb);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nb),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        boost::container::small_vector<int, 32> idx;
        for (std::size_t b = r.begin(); b != r.end(); ++b)
        {
          std::vector<int>& u = unions[b];
          for (const Element_type& e : parts[b])
          {
            zone_grid_indices(e, op, c3t3, idx);
            u.insert(u.end(), idx.begin(), idx.end());
          }
          std::sort(u.begin(), u.end());
          u.erase(std::unique(u.begin(), u.end()), u.end());
        }
      });

    // Cap: a bucket claiming a large share of everything the phase touches
    // would serialize the phase, so it is not even attempted.
    std::size_t total = 0;
    for (const auto& u : unions) total += u.size();
    const std::size_t cap = (nb == 0) ? 0 : (2 * total) / (nb ? nb : 1);

    const auto& tr = c3t3.triangulation();
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nb, 1),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t b = r.begin(); b != r.end(); ++b)
        {
          const std::vector<int>& u = unions[b];
          bool whole = !u.empty() && u.size() <= cap;
          if (whole)
          {
            for (const int gi : u)              // ascending: no cyclic wait
              if (!tr.try_lock_grid_index(gi)) { whole = false; break; }
          }
          if (whole)
          {
            for (const Element_type& e : parts[b])
              op.execute_operation(e, c3t3);
            tr.unlock_all_elements();
          }
          else
          {
            tr.unlock_all_elements();
            for (const Element_type& e : parts[b])
              apply_one(e, op, c3t3);
          }
        }
      });
  }

  /**
  * B4 -- a 27-colour sweep with no per-element locking.
  *
  * Colour each element by the parity-mod-3 of its lock-grid cell. Two elements
  * of the same colour sit in grid cells whose indices differ by a multiple of
  * three in every axis, hence by at least three in some axis whenever they
  * differ at all -- so the 27-neighbourhoods of those cells are disjoint, and
  * two zones confined to their own neighbourhoods cannot overlap. 27 rather
  * than the 8 of a red/black scheme: at distance two the neighbourhoods still
  * share a plane of cells, so eight colours would NOT be sound here.
  *
  * An element whose zone escapes its cell's 27-neighbourhood goes to a
  * remainder set run afterwards with ordinary locking, so soundness does not
  * depend on the containment being universal -- only on it being tested.
  */
  static void run_coloured(std::vector<std::vector<Element_type> >& parts,
                           Operation& op, C3t3& c3t3)
  {
    const auto& tr = c3t3.triangulation();
    std::vector<std::vector<Element_type> > colours(27);
    std::vector<Element_type> remainder;

    boost::container::small_vector<int, 32> idx;
    for (const auto& part : parts)
      for (const Element_type& e : part)
      {
        boost::container::small_vector<Vertex_handle, 2> lv;
        op.elision_vertices(e, lv);
        if (lv.empty()) { remainder.push_back(e); continue; }
        const std::array<int, 3> c = tr.lock_grid_indices3(lv[0]->point());
        if (c[0] < 0) { remainder.push_back(e); continue; }

        zone_grid_indices(e, op, c3t3, idx);
        bool contained = !idx.empty();
        for (const int gi : idx)
        {
          const std::array<int, 3> g = unflatten(tr, gi);
          if (std::abs(g[0] - c[0]) > 1 || std::abs(g[1] - c[1]) > 1
           || std::abs(g[2] - c[2]) > 1) { contained = false; break; }
        }
        if (!contained) { remainder.push_back(e); continue; }
        colours[(c[0] % 3) + 3 * (c[1] % 3) + 9 * (c[2] % 3)].push_back(e);
      }

    // Within one colour the grid cells are pairwise non-adjacent, so zones in
    // DIFFERENT cells cannot overlap -- but two elements in the SAME cell can,
    // and they must therefore run in the same task. Grouping by cell is what
    // gives that; partitioning the colour by kd median does not, and running
    // it that way segfaulted.
    for (std::vector<Element_type>& colour : colours)
    {
      if (colour.empty())
        continue;
      std::unordered_map<int, std::vector<Element_type> > by_cell;
      for (const Element_type& e : colour)
      {
        boost::container::small_vector<Vertex_handle, 2> lv;
        op.elision_vertices(e, lv);
        by_cell[tr.lock_grid_index(lv[0]->point())].push_back(e);
      }
      std::vector<std::vector<Element_type> > groups;
      groups.reserve(by_cell.size());
      for (auto& kv : by_cell)
        groups.push_back(std::move(kv.second));

      tbb::parallel_for_each(groups.begin(), groups.end(),
        [&](const std::vector<Element_type>& group)
        {
          for (const Element_type& e : group)
            op.execute_operation(e, c3t3);
          c3t3.triangulation().unlock_all_elements();
        });
    }

    if (!remainder.empty())
    {
      std::vector<std::vector<Element_type> > sub = kd_partition(remainder, op);
      const Interior_flags none(sub.size()); // one empty flag row per bucket
      run_flagged(sub, none, op, c3t3);
    }
  }

  static std::array<int, 3> unflatten(const typename C3t3::Triangulation& tr,
                                      int gi)
  {
    // grid_index = z*n^2 + y*n + x; recovered through a probe of the axis count
    const int n = grid_axis(tr);
    const int z = gi / (n * n);
    gi -= z * n * n;
    const int y = gi / n;
    return { gi - y * n, y, z };
  }

  static int grid_axis(const typename C3t3::Triangulation& tr)
  {
    return tr.lock_grid_num_cells_per_axis();
  }

  static void run_parts(std::vector<std::vector<Element_type> >& parts,
                        Operation& op, C3t3& c3t3)
  {
    // CGAL_TR_ELISION_MODE selects one of the elision strategies. They are
    // alternatives, not layers, so a single integer makes it impossible to
    // combine two by accident. 0 off, 1 R19 as measured, 2 B1, 3 B2, 4 B3,
    // 5 B4.
    switch (Parallel_tuning::get().elision_mode)
    {
    case 2: return run_flagged(parts, classify_flags_vertex(parts, op, c3t3), op, c3t3);
    case 3: return run_flagged(parts, classify_flags_grid(parts, op, c3t3), op, c3t3);
    case 4: return run_bucket_scoped(parts, op, c3t3);
    case 5: return run_coloured(parts, op, c3t3);
    default: break;
    }

    const Interior_set interior = classify_interior(parts, op, c3t3);
    if (elision_stats())
    {
      std::size_t n_int = 0, n_tot = 0;
      for (const std::vector<Element_type>& part : parts)
        for (const Element_type& e : part)
        { ++n_tot; if (is_interior(e, op, interior)) ++n_int; }
      std::fprintf(stderr, "[elision] %-34s interior %8zu | locked %8zu | %5.1f%% freed\n",
                   op.operation_name().c_str(), n_int, n_tot - n_int,
                   n_tot ? 100.0 * double(n_int) / double(n_tot) : 0.0);
    }
#ifdef CGAL_TR_ELISION_CHECK
    check_interior_classification(parts, op, c3t3, interior);
#endif

    const bool defer = Parallel_tuning::get().defer_on_conflict;

    tbb::parallel_for_each(parts.begin(), parts.end(),
                           [&](const std::vector<Element_type>& part)
                           {
                             std::vector<Element_type> deferred;
                             for (const Element_type& element : part)
                             {
                               if (is_interior(element, op, interior))
                                 apply_one_unlocked(element, op, c3t3);
                               else if (defer)
                                 apply_one_or_defer(element, op, c3t3, deferred);
                               else
                                 apply_one(element, op, c3t3);
                             }
                             // backstop: whatever was contended twice is run
                             // now, with the unbounded loop
                             for (const Element_type& element : deferred)
                               apply_one(element, op, c3t3);
                           });
  }

#ifdef CGAL_TR_ELISION_CHECK
  /**
  * Checks the property the elision rests on, independently of how it was
  * derived. The lock probe cannot do this: for an elided element "no lock
  * held" is correct by design, so a misclassified element and a correctly
  * elided one look identical to it.
  *
  * For every element classified interior, this walks the ACTUAL zone -- the
  * cells incident to each of its locked vertices -- and requires that every
  * vertex of every one of those cells is either untagged (owns no element, so
  * locks nothing) or owned by this same bucket. If that fails, some other
  * bucket can reach the zone and running it unlocked is a race.
  *
  * Serial and timing independent: it runs before the parallel phase, on a mesh
  * nobody is modifying, so it fires on the first run that misclassifies rather
  * than on the first run that collides.
  */
  static void check_interior_classification(
      const std::vector<std::vector<Element_type> >& parts,
      const Operation& op, const C3t3& c3t3, const Interior_set& interior)
  {
    if (interior.empty())
      return;

    // which bucket each element's locked vertices belong to
    std::unordered_map<Vertex_handle, std::size_t> bucket_of;
    boost::container::small_vector<Vertex_handle, 2> lv;
    for (std::size_t b = 0; b < parts.size(); ++b)
      for (const Element_type& e : parts[b])
      {
        lv.clear();
        op.locked_vertices(e, lv);
        for (const Vertex_handle& v : lv)
          bucket_of.emplace(v, b);
      }

    std::size_t checked = 0, violations = 0;
    const auto& tr = c3t3.triangulation();
    for (std::size_t b = 0; b < parts.size(); ++b)
      for (const Element_type& e : parts[b])
      {
        if (!is_interior(e, op, interior))
          continue;
        ++checked;
        lv.clear();
        op.locked_vertices(e, lv);
        for (const Vertex_handle& v : lv)
        {
          std::vector<typename C3t3::Triangulation::Cell_handle> star;
          tr.incident_cells(v, std::back_inserter(star));
          for (const auto& c : star)
            for (int k = 0; k < 4; ++k)
            {
              const auto it = bucket_of.find(c->vertex(k));
              if (it != bucket_of.end() && it->second != b)
              {
                ++violations;
                std::cerr << "[elision] VIOLATION: element in bucket " << b
                          << " classified interior, but a cell of its zone has"
                             " a vertex owned by bucket " << it->second
                          << std::endl;
              }
            }
        }
      }
    std::size_t total = 0;
    for (const auto& part : parts) total += part.size();
    std::cerr << "[elision] " << checked << " of " << total
              << " elements interior ("
              << (total ? 100.0 * double(checked) / double(total) : 0.0)
              << "%), " << violations << " zone violations, "
              << parts.size() << " buckets" << std::endl;
  }
#endif

  static void run_buckets(Buckets& buckets, Operation& op, C3t3& c3t3)
  {
    tbb::parallel_for_each(buckets.begin(), buckets.end(),
                           [&](const std::pair<const Grid_cell_index,
                                               std::vector<Element_type>>& bucket)
                           {
                             for (const Element_type& element : bucket.second)
                               apply_one(element, op, c3t3);
                           });
  }

  /**
  * Partitions the elements into equal-count buckets and runs one per task.
  * Elements in a bucket are close together, so a thread that has just locked
  * one zone tends to find the next adjacent rather than contended, and two
  * threads are usually working in different regions. Buckets hold equal
  * COUNTS rather than equal volumes, which is what keeps the threads busy for
  * the same length of time: grouping by position instead was measured and
  * superseded, because an unevenly dense mesh gave one oversized bucket that
  * serialized the phase.
  */
  /** D3 -- Hilbert order, then equal-size contiguous chunks. */
  static std::vector<std::vector<Element_type> >
  hilbert_partition(std::vector<Element_type>& candidates, const Operation& op)
  {
    using GT = typename C3t3::Triangulation::Geom_traits;
    using Point_3 = typename GT::Point_3;

    const std::size_t n = candidates.size();
    std::vector<Point_3> pts;
    pts.reserve(n);
    for (const Element_type& e : candidates)
      pts.push_back(op.point_on_element(e));

    std::vector<std::ptrdiff_t> ids(n);
    std::iota(ids.begin(), ids.end(), std::ptrdiff_t(0));

    // Equal counts by construction, as the shipped equal-count kd partition
    // gives -- but the chunks are also COMPACT. A kd tree with 16 leaves cuts
    // four times along axis-aligned directions and its boxes are elongated
    // whenever the per-axis split counts differ; elongated regions have large
    // surface area, hence a large interface with their neighbours, and the
    // interface is what produces lock conflicts and blocks elision. Hilbert
    // chunks are compact by the curve's locality property, and elements
    // adjacent in the order are adjacent in space, so a thread walks the mesh
    // coherently inside a bucket. Cost is one sort against one nth_element
    // recursion, so this is close to cost-neutral -- unlike D1.
    CGAL::Spatial_sort_traits_adapter_3<GT, const Point_3*> traits(pts.data());
    CGAL::hilbert_sort(ids.begin(), ids.end(), traits);

    const std::size_t target = bucket_target();
    std::vector<std::vector<Element_type> > parts;
    parts.reserve(target);
    for (std::size_t b = 0; b < target; ++b)
    {
      const std::size_t lo = (n * b) / target;
      const std::size_t hi = (n * (b + 1)) / target;
      if (hi <= lo)
        continue;
      std::vector<Element_type> part;
      part.reserve(hi - lo);
      for (std::size_t i = lo; i < hi; ++i)
        part.push_back(candidates[static_cast<std::size_t>(ids[i])]);
      parts.push_back(std::move(part));
    }
    return parts;
  }

  static std::size_t bucket_target()
  {
    return static_cast<std::size_t>(Parallel_tuning::get().buckets_per_thread)
         * static_cast<std::size_t>((std::max)(1, tbb::this_task_arena::max_concurrency()));
  }

  /**
  * D1 -- METIS on the candidate adjacency graph.
  *
  * Nodes are the candidates; two candidates are joined when their lock zones
  * can overlap, which is detected here by them sharing a lock-grid cell.
  * METIS_PartGraphKway then minimises the number of edges cut subject to a
  * balance constraint -- that is, it minimises exactly the INTERFACE between
  * buckets, which the equal-count kd partition does not optimise at all. Fewer
  * cross-bucket adjacencies means fewer lock conflicts and a larger elidable
  * interior, so this is also the multiplier on Group B.
  *
  * The cost is the whole risk: building the graph is Theta(candidates x zone)
  * and METIS is not free, and this runs per operation per phase. Returns false
  * if METIS is unavailable or fails, and the caller falls back.
  */
  static bool metis_partition(std::vector<Element_type>& candidates,
                              const Operation& op, const C3t3& c3t3,
                              std::vector<std::vector<Element_type> >& parts)
  {
    const Metis_api& api = Metis_api::get();
    const std::size_t n = candidates.size();
    const std::size_t target = bucket_target();
    if (!api.ok || n < target * 2 || target < 2)
      return false;

    // candidates that share a lock-grid cell are adjacent
    std::unordered_map<int, std::vector<int> > by_cell;
    std::vector<boost::container::small_vector<int, 16> > zones(n);
    const auto& tr = c3t3.triangulation();
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, n),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t i = r.begin(); i != r.end(); ++i)
          zone_grid_indices(candidates[i], op, c3t3, zones[i]);
      });
    for (std::size_t i = 0; i < n; ++i)
    {
      if (zones[i].empty())
      {
        const int gi = tr.lock_grid_index(op.point_on_element(candidates[i]));
        if (gi >= 0)
          zones[i].push_back(gi);
      }
      for (const int gi : zones[i])
        by_cell[gi].push_back(static_cast<int>(i));
    }

    std::vector<std::unordered_set<int> > adj(n);
    for (const auto& kv : by_cell)
    {
      const std::vector<int>& members = kv.second;
      if (members.size() > 64)   // a cell everybody touches carries no signal
        continue;
      for (std::size_t a = 0; a < members.size(); ++a)
        for (std::size_t b = a + 1; b < members.size(); ++b)
        {
          adj[members[a]].insert(members[b]);
          adj[members[b]].insert(members[a]);
        }
    }

    std::vector<int> xadj(n + 1, 0), adjncy;
    adjncy.reserve(n * 8);
    for (std::size_t i = 0; i < n; ++i)
    {
      for (const int j : adj[i])
        adjncy.push_back(j);
      xadj[i + 1] = static_cast<int>(adjncy.size());
    }
    if (adjncy.empty())
      return false;

    int nvtxs = static_cast<int>(n), ncon = 1, nparts = static_cast<int>(target);
    int edgecut = 0;
    std::vector<int> part(n, 0), options(40, -1);
    api.set_opts(options.data());
    const int rc = api.part_kway(&nvtxs, &ncon, xadj.data(), adjncy.data(),
                                 nullptr, nullptr, nullptr, &nparts,
                                 nullptr, nullptr, options.data(),
                                 &edgecut, part.data());
    if (rc != 1)   // METIS_OK
      return false;
    for (const int p : part)
      if (p < 0 || p >= nparts)
        return false;   // not a partition we can trust

    parts.assign(static_cast<std::size_t>(nparts), std::vector<Element_type>());
    for (std::size_t i = 0; i < n; ++i)
      parts[static_cast<std::size_t>(part[i])].push_back(candidates[i]);
    parts.erase(std::remove_if(parts.begin(), parts.end(),
                  [](const std::vector<Element_type>& p) { return p.empty(); }),
                parts.end());
    return parts.size() > 1;
  }

  /** D2 -- build the shared subdivision's split planes, or reuse them. */
  static void build_shared_partition(std::vector<Element_type>& candidates,
                                     const Operation& op)
  {
    Shared_kd_partition& sp = Shared_kd_partition::get();
    std::size_t leaves = 1;
    while (leaves < bucket_target())
      leaves *= 2;                      // perfect tree: a power of two
    sp.nodes.assign(leaves - 1, Shared_kd_partition::Node{});
    sp.leaves = leaves;

    std::vector<Element_type> scratch(candidates);
    record_split(scratch.begin(), scratch.end(), op, 0, leaves, sp);
    sp.valid = true;
  }

  static void record_split(typename std::vector<Element_type>::iterator first,
                           typename std::vector<Element_type>::iterator last,
                           const Operation& op, std::size_t node,
                           std::size_t span, Shared_kd_partition& sp)
  {
    if (span <= 1 || node >= sp.nodes.size())
      return;
    const std::size_t n = static_cast<std::size_t>(std::distance(first, last));
    if (n < 2)
    {
      record_split(first, last, op, 2 * node + 1, span / 2, sp);
      record_split(first, last, op, 2 * node + 2, span / 2, sp);
      return;
    }
    double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
    double hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
    for (auto it = first; it != last; ++it)
    {
      const auto p = op.point_on_element(*it);
      const double c[3] = { CGAL::to_double(p.x()), CGAL::to_double(p.y()),
                            CGAL::to_double(p.z()) };
      for (int a = 0; a < 3; ++a)
      { lo[a] = (std::min)(lo[a], c[a]); hi[a] = (std::max)(hi[a], c[a]); }
    }
    int axis = 0;
    for (int a = 1; a < 3; ++a)
      if (hi[a] - lo[a] > hi[axis] - lo[axis]) axis = a;

    const auto coord = [&](const Element_type& e)
    {
      const auto p = op.point_on_element(e);
      return (axis == 0) ? CGAL::to_double(p.x())
           : (axis == 1) ? CGAL::to_double(p.y()) : CGAL::to_double(p.z());
    };
    const auto mid = first + static_cast<std::ptrdiff_t>(n / 2);
    std::nth_element(first, mid, last,
      [&](const Element_type& a, const Element_type& b) { return coord(a) < coord(b); });

    sp.nodes[node].axis = axis;
    sp.nodes[node].value = coord(*mid);
    record_split(first, mid, op, 2 * node + 1, span / 2, sp);
    record_split(mid, last, op, 2 * node + 2, span / 2, sp);
  }

  static std::vector<std::vector<Element_type> >
  shared_partition(std::vector<Element_type>& candidates, const Operation& op)
  {
    Shared_kd_partition& sp = Shared_kd_partition::get();
    if (!sp.valid)
      build_shared_partition(candidates, op);
    std::vector<std::vector<Element_type> > parts(sp.leaves);
    for (const Element_type& e : candidates)
    {
      const auto p = op.point_on_element(e);
      parts[sp.leaf_of(CGAL::to_double(p.x()), CGAL::to_double(p.y()),
                       CGAL::to_double(p.z()))].push_back(e);
    }
    parts.erase(std::remove_if(parts.begin(), parts.end(),
                  [](const std::vector<Element_type>& p) { return p.empty(); }),
                parts.end());
    return parts;
  }

  static void run_unordered(std::vector<Element_type>& candidates,
                            Operation& op, C3t3& c3t3)
  {
    const Parallel_tuning& t = Parallel_tuning::get();
    std::vector<std::vector<Element_type> > parts;

    if (t.partitioner == 2 && metis_partition(candidates, op, c3t3, parts))
      ;                                            // D1
    else if (t.partitioner == 1)
      parts = hilbert_partition(candidates, op);   // D3
    else if (t.reuse_partition)
      parts = shared_partition(candidates, op);    // D2
    else
      parts = kd_partition(candidates, op);

    // D4c -- largest bucket first (LPT). Equal-count buckets hold equal
    // COUNTS, but the cost per element is not equal (a collapse that fails
    // costs far less than one that succeeds), so bucket runtimes still vary
    // and the phase waits for the slowest. Submitting the big ones first is
    // the classic remedy and costs one sort of ~16 pointers.
    if (t.bucket_schedule == 1)
      std::stable_sort(parts.begin(), parts.end(),
                       [](const std::vector<Element_type>& a,
                          const std::vector<Element_type>& b)
                       { return a.size() > b.size(); });

    run_parts(parts, op, c3t3);
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
