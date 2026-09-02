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

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/blocked_range.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/task_arena.h>
#endif

#include <algorithm>
#include <atomic>
#include <iterator>
#include <string>
#include <thread>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <random>
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
    while (!op.lock_zone(element, c3t3))
    {
      c3t3.triangulation().unlock_all_elements();
      std::this_thread::yield();
    }
    op.execute_operation(element, c3t3);
    c3t3.triangulation().unlock_all_elements();
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
  /**
  * Both arms live in one binary, selected at run time, so that comparing them
  * compares the grouping and not two compilations (POLICY 0.2). Unset means
  * grouped; `CGAL_TR_BUCKET_UNORDERED=0` restores the shuffle.
  */
  static bool bucket_unordered_enabled()
  {
    static const bool enabled = []
      {
        const char* const e = std::getenv("CGAL_TR_BUCKET_UNORDERED");
        return (e == nullptr) || (std::atoi(e) != 0);
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

  // The shuffle the grouping replaced: spread the threads over the
  // triangulation by chance rather than by construction.
  static void run_unordered_shuffled(std::vector<Element_type>& candidates,
                                     Operation& op, C3t3& c3t3)
  {
    std::mt19937 gen(std::random_device{}());
    std::shuffle(candidates.begin(), candidates.end(), gen);

    tbb::parallel_for_each(candidates,
                           [&](const Element_type& element)
                           {
                             apply_one(element, op, c3t3);
                           });
  }

  /**
  * Groups the elements by grid cell and runs one bucket per task. Elements in
  * a bucket are close together, so a thread that has just locked one zone is
  * likely to find the next adjacent rather than contended, and two threads are
  * usually working in different regions.
  */
  static void run_unordered(std::vector<Element_type>& candidates,
                            Operation& op, C3t3& c3t3)
  {
    if (!bucket_unordered_enabled())
      return run_unordered_shuffled(candidates, op, c3t3);

    Buckets buckets = bucket_by_grid(candidates, op, c3t3);
    if (buckets.empty()) // degenerate bounding box
      return run_unordered_shuffled(candidates, op, c3t3);

    run_buckets(buckets, op, c3t3);
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
