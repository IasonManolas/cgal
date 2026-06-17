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
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Iterator_range.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/task_group.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_vector.h>
#include <tbb/task_arena.h>
#include <functional>
#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/functional/hash.hpp>
#endif

#ifndef LOCK_GRID_SIZE
#  define LOCK_GRID_SIZE 128 // grid cells per axis for the spatial lock grid
#endif

#ifndef LB_BUCKETS_PER_THREAD
// Unordered ops partition the candidate set into ~LB_BUCKETS_PER_THREAD * nthreads
// equal-count, spatially-compact kd-buckets so TBB work-stealing balances load.
#  define LB_BUCKETS_PER_THREAD 4
#endif

#include <atomic>
#include <algorithm>
#include <random>
#include <vector>
#include <boost/container/small_vector.hpp>
#include <thread>
#include <string>
#include <utility>
#include <iostream>
#include <fstream>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

template <typename C3t3_, typename ElementType_, typename ElementSource_, typename LockElementType_>
class ElementaryOperation
{
public:
  using C3t3 = C3t3_;
  using Triangulation = typename C3t3::Triangulation;
  using ElementType = ElementType_;
  using ElementSource = ElementSource_;
  using LockElementType = LockElementType_;
  using Lock_zone = std::vector<LockElementType>;

  ElementaryOperation() = default;
  virtual ~ElementaryOperation() = default;

  // Compile-time processing strategy trait
  static constexpr bool ordered_processing = false;
  virtual bool requires_ordered_processing() const { return false; }
  // Pure element logic methods
  virtual ElementSource get_element_source(const C3t3& c3t3) const = 0;
  virtual bool lock_zone(const ElementType& e, const C3t3& c3t3) const = 0;
  virtual bool execute_operation(const ElementType& e, C3t3& c3t3) = 0;
  virtual std::string operation_name() const = 0;

  // The vertices whose incident-cell 1-rings lock_zone() locks. Used by
  // apply_unordered_processing for interior/boundary lock elision: an element all of
  // whose locked vertices are bucket-interior (never shared with another
  // concurrently-processed bucket) can run lock-free. Default: empty -> the element
  // is conservatively treated as BOUNDARY (always locked). Flip and smooth override.
  using Vertex_handle = typename Triangulation::Vertex_handle;
  virtual void locked_vertices(const ElementType& e,
                               boost::container::small_vector<Vertex_handle, 2>& out) const
  { (void)e; (void)out; }
};

// Base class for operation execution strategies
template <typename Operation> class ElementaryOperationExecution
{
public:
  using C3t3 = typename Operation::C3t3;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;
  using ElementType = typename Operation::ElementType;

  ElementaryOperationExecution() = default;
  ElementaryOperationExecution(const ElementaryOperationExecution&) = default;
  ElementaryOperationExecution(ElementaryOperationExecution&&) = default;
  ElementaryOperationExecution& operator=(const ElementaryOperationExecution&) = default;
  ElementaryOperationExecution& operator=(ElementaryOperationExecution&&) = default;
  virtual ~ElementaryOperationExecution() = default;

  virtual std::vector<ElementType> collect_candidates(const Operation& op, const C3t3& c3t3) const = 0;

  virtual bool apply_operation_on_elements(std::vector<ElementType>& elements, Operation& op, C3t3& c3t3) = 0;

  // Main execution method that orchestrates the operation
  virtual bool execute(Operation& op, C3t3& c3t3) {
    // Collect candidate elements
    auto candidates = collect_candidates(op, c3t3);

    if(candidates.empty()) {
      return false;
    }

#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
    // Ensure lock data structure is initialized for parallel mode
    using Tr = typename C3t3::Triangulation;
    using Concurrency = typename Tr::Concurrency_tag;
    if constexpr(std::is_convertible<Concurrency, CGAL::Parallel_tag>::value)
      ensure_lock_data_structure_initialized(c3t3);
#endif

    // Apply operations on elements
    return apply_operation_on_elements(candidates, op, c3t3);
  }

private:
  // Ensure the lock data structure is initialized when needed
  void ensure_lock_data_structure_initialized(C3t3& c3t3) {
    auto& triangulation = c3t3.triangulation();
    if(!triangulation.get_lock_data_structure()) {
      // Create a lock data structure using the C3T3's bounding box

#ifndef USE_TAG_NON_BLOCKING
      static typename C3t3::Triangulation::Lock_data_structure lock_ds(c3t3.bbox(), // Use C3T3's bounding box
                                                                       LOCK_GRID_SIZE // grid size (overridable)
      );
#else
      static CGAL::Spatial_lock_grid_3<Tag_non_blocking> lock_ds(c3t3.bbox(), LOCK_GRID_SIZE);
#endif
      triangulation.set_lock_data_structure(&lock_ds);
    }
  }
};

// Sequential execution implementation
template <typename Operation>
class ElementaryOperationExecutionSequential : public ElementaryOperationExecution<Operation>
{
public:
  using Base = ElementaryOperationExecution<Operation>;
  using C3t3 = typename Operation::C3t3;
  using ElementType = typename Operation::ElementType;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;

public:
  ElementaryOperationExecutionSequential() = default;
  ElementaryOperationExecutionSequential(const ElementaryOperationExecutionSequential&) = default;
  ElementaryOperationExecutionSequential(ElementaryOperationExecutionSequential&&) = default;
  ElementaryOperationExecutionSequential& operator=(const ElementaryOperationExecutionSequential&) = default;
  ElementaryOperationExecutionSequential& operator=(ElementaryOperationExecutionSequential&&) = default;

  std::vector<ElementType> collect_candidates(const Operation& op, const C3t3& c3t3) const override {
    auto elements = op.get_element_source(c3t3);
    std::vector<ElementType> candidates;
    // candidates.reserve(elements.size());
    // std::copy(elements.begin(), elements.end(), std::back_inserter(candidates));

    if constexpr(std::is_same<decltype(elements), std::vector<ElementType>>::value) {
      candidates = std::move(elements); // directly move if vector
    } else {
      candidates.reserve(std::distance(elements.begin(), elements.end()));
      std::copy(elements.begin(), elements.end(), std::back_inserter(candidates));
    }
    return candidates;
  }

  bool apply_operation_on_elements(std::vector<ElementType>& elements, Operation& op, C3t3& c3t3) override {
#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
    std::size_t num_successful_locks = 0;
#endif

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Executing operation sequentially: " << op.operation_name() << "...";
#endif

    for(const auto& element : elements) {
      if(op.execute_operation(element, c3t3)) {
#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
        num_successful_locks++;
#endif
      }
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << " done num_elements:" << elements.size()
#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
              << ", num_successful_locks : " << num_successful_locks
#endif
              << std::endl;
#endif

    return true;
  }
};

template <typename Operation>
class ElementaryOperationExecutionParallel : public ElementaryOperationExecution<Operation>
{
public:
  using Base = ElementaryOperationExecution<Operation>;
  using C3t3 = typename Operation::C3t3;
  using ElementType = typename Operation::ElementType;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;

private:
  size_t m_batch_size;

public:
  ElementaryOperationExecutionParallel(size_t batch_size = 64)
      : m_batch_size(batch_size) {}
  ElementaryOperationExecutionParallel(const ElementaryOperationExecutionParallel&) = default;
  ElementaryOperationExecutionParallel(ElementaryOperationExecutionParallel&&) = default;
  ElementaryOperationExecutionParallel& operator=(const ElementaryOperationExecutionParallel&) = default;
  ElementaryOperationExecutionParallel& operator=(ElementaryOperationExecutionParallel&&) = default;

  std::vector<ElementType> collect_candidates(const Operation& op, const C3t3& c3t3) const override {
    auto elements = op.get_element_source(c3t3);
    std::vector<ElementType> candidates;
    if constexpr(std::is_same<decltype(elements), std::vector<ElementType>>::value) {
      candidates = std::move(elements); // directly move if vector
    } else {
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
      tbb::combinable<std::vector<ElementType>> local_candidates;
      tbb::parallel_for_each(elements.begin(), elements.end(),
                             [&](const ElementType& element) { local_candidates.local().push_back(element); });
      local_candidates.combine_each([&](const std::vector<ElementType>& local) {
        candidates.insert(candidates.end(), local.begin(), local.end());
      });
#endif
    }
    return candidates;
  }

private:
  // Ordered processing using spatial bucketing (for operations like EdgeSplit, EdgeCollapse)
  bool apply_ordered_processing(
    std::vector<ElementType>& elements,
    Operation& op,
    C3t3& c3t3,
    const CGAL::Bbox_3& bb)
  {
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB

#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
    std::atomic<size_t> num_successful_locks = 0;
    std::atomic<size_t> num_failed_locks = 0;
#endif

    const double min_sq_dim = (std::min)(CGAL::square(bb.xmax() - bb.xmin()),
                                (std::min)(CGAL::square(bb.ymax() - bb.ymin()),
                                           CGAL::square(bb.zmax() - bb.zmin())));
    const float grid_cell_size = 0.5f * CGAL::approximate_sqrt(min_sq_dim);
    const float inv_cell_size = 1.f / grid_cell_size;

    struct Grid_cell_index
    {
      int i, j, k;
      bool operator==(const Grid_cell_index& other) const
      { return i == other.i && j == other.j && k == other.k; }
    };

    struct Grid_cell_index_hasher
    {
      std::size_t operator()(const Grid_cell_index& ci) const
      {
        return ((std::hash<int>()(ci.i) ^ (std::hash<int>()(ci.j) << 1)) >> 1)
               ^ (std::hash<int>()(ci.k) << 1);
      }
    };

    auto compute_cell_index
      = [&](const typename C3t3::Triangulation::Geom_traits::Point_3& p) -> Grid_cell_index
        {
          return {static_cast<int>(std::floor(p.x() * inv_cell_size)),
                  static_cast<int>(std::floor(p.y() * inv_cell_size)),
                  static_cast<int>(std::floor(p.z() * inv_cell_size))};
        };

    std::unordered_map<Grid_cell_index, std::vector<ElementType>, Grid_cell_index_hasher> spatialBuckets;
    for(const auto& e : elements)
    {
      Grid_cell_index idx = compute_cell_index(op.point_on_element(e));
      spatialBuckets[idx].push_back(e);
    }

    tbb::parallel_for_each(spatialBuckets.begin(), spatialBuckets.end(),
      [&](const std::pair<const Grid_cell_index, std::vector<ElementType>>& bucket)
      {
        for(const auto& element : bucket.second)
        {
          while(!op.lock_zone(element, c3t3))
          {
            c3t3.triangulation().unlock_all_elements();
            std::this_thread::yield();
#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
            num_failed_locks++;
#endif
          }
          op.execute_operation(element, c3t3);
#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
          num_successful_locks++;
#endif
          c3t3.triangulation().unlock_all_elements();
        }
      });

#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
    {
      const char* csv_path = "lock_stats.csv";
      bool need_header = true;
      {
        std::ifstream ifs(csv_path);
        need_header = !ifs.good() || (ifs.peek() == std::ifstream::traits_type::eof());
      }
      std::ofstream ofs(csv_path, std::ios::app);
      if(need_header)
        ofs << "operation,num_successful_locks,num_failed_locks,lock_success_rate(%)" << std::endl;
      const size_t attempts = static_cast<size_t>(num_successful_locks + num_failed_locks);
      const double lock_success_rate =
          (attempts == 0) ? 0.0 : 100.0 * static_cast<double>(num_successful_locks) / static_cast<double>(attempts);
      ofs << op.operation_name() << "," << num_successful_locks << "," << num_failed_locks
          << "," << lock_success_rate << std::endl;
    }
#endif
#endif // concurrent

    return true;
  }

  // Unordered processing (VertexSmooth, EdgeFlip): no ordering constraint. Partition
  // the candidate set into equal-count, spatially-compact buckets via recursive
  // median (kd-tree) splitting, then process buckets concurrently with each bucket
  // sequential on one thread. Equal counts -> balanced load (no idle tail); compact
  // boxes -> lock conflicts stay confined to bucket boundaries (interiors never
  // contend). Replaces the previous uniform 0.5*sqrt(min bbox dim) grid, whose few
  // lopsided buckets starved workers at the tail.
  bool apply_unordered_processing(
    std::vector<ElementType>& elements,
    Operation& op,
    C3t3& c3t3,
    const CGAL::Bbox_3& bb)
  {
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB

#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
    std::atomic<size_t> num_successful_locks = 0;
    std::atomic<size_t> num_failed_locks = 0;
    std::atomic<size_t> num_elided_locks = 0;
#endif

    (void)bb; // partition is data-driven (kd-tree over element points), not bbox-grid

    // Precompute element points once (double coords) for the median split,
    // in parallel: point_on_element + to_double over the whole candidate set.
    struct EP { double c[3]; ElementType e; };
    std::vector<EP> eps(elements.size());
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, elements.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for(std::size_t i = r.begin(); i != r.end(); ++i)
        {
          const auto p = op.point_on_element(elements[i]);
          eps[i] = EP{{CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())}, elements[i]};
        }
      });

    // Target bucket size: enough buckets per thread for work-stealing slack.
    const std::size_t nthreads =
      static_cast<std::size_t>((std::max)(1, tbb::this_task_arena::max_concurrency()));
    const std::size_t target =
      (std::max)(std::size_t(1), eps.size() / (std::size_t(LB_BUCKETS_PER_THREAD) * nthreads));

    // Recursive median split -> equal-count index ranges [lo,hi). The two halves
    // of each split are independent, so above a threshold they recurse as
    // concurrent tasks (task_group); this keeps the build off the serial critical
    // path (serial chain drops from O(N log B) to ~O(2N)). Below the threshold the
    // recursion runs serially to avoid task overhead.
    constexpr std::size_t PAR_SPLIT_THRESHOLD = 8192;
    tbb::concurrent_vector<std::pair<std::size_t, std::size_t>> buckets;
    std::function<void(std::size_t, std::size_t)> split_range =
      [&](std::size_t lo, std::size_t hi)
      {
        if(hi - lo <= target)
        {
          if(hi > lo)
            buckets.emplace_back(lo, hi);
          return;
        }
        // Split along the axis of largest spatial extent at the median element.
        double mn[3] = {eps[lo].c[0], eps[lo].c[1], eps[lo].c[2]};
        double mx[3] = {eps[lo].c[0], eps[lo].c[1], eps[lo].c[2]};
        for(std::size_t i = lo + 1; i < hi; ++i)
          for(int a = 0; a < 3; ++a)
          {
            mn[a] = (std::min)(mn[a], eps[i].c[a]);
            mx[a] = (std::max)(mx[a], eps[i].c[a]);
          }
        int axis = 0;
        double best = mx[0] - mn[0];
        for(int a = 1; a < 3; ++a)
          if(mx[a] - mn[a] > best) { best = mx[a] - mn[a]; axis = a; }

        const std::size_t mid = lo + (hi - lo) / 2;
        std::nth_element(eps.begin() + lo, eps.begin() + mid, eps.begin() + hi,
                         [axis](const EP& x, const EP& y) { return x.c[axis] < y.c[axis]; });
        if(hi - lo > PAR_SPLIT_THRESHOLD)
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

    // --- interior/boundary classification for lock elision (topological) ---
    // Tag each element-endpoint vertex with its bucket id (MIXED if >=2 buckets own
    // it); then flag every vertex incident to a cell whose tagged vertices span >=2
    // buckets ("boundary-touching"). An element is INTERIOR iff none of its locked
    // vertices is boundary-touching => its whole 1-ring lock zone is single-bucket, so
    // no other concurrently-processed bucket can touch it and it may run lock-free.
    // flip/smooth create/destroy no vertices, so this classification stays valid for
    // the whole phase. Empty locked_vertices() => the element falls through to the
    // locked path (the safe default).
    using Vh = typename C3t3::Triangulation::Vertex_handle;
    const std::size_t MIXED_BUCKET = static_cast<std::size_t>(-1);
    boost::concurrent_flat_map<Vh, std::size_t, boost::hash<Vh>> vertex_bucket;
    vertex_bucket.reserve(eps.size() * 2);

    // Pass 1: scatter element-endpoint vertex -> bucket id (MIXED on conflict).
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, buckets.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        boost::container::small_vector<Vh, 2> lv;
        for(std::size_t bi = r.begin(); bi != r.end(); ++bi)
          for(std::size_t i = buckets[bi].first; i < buckets[bi].second; ++i)
          {
            lv.clear();
            op.locked_vertices(eps[i].e, lv);
            for(const Vh& v : lv)
              vertex_bucket.insert_or_visit(std::make_pair(v, bi),
                [bi, MIXED_BUCKET](std::pair<const Vh, std::size_t>& kv)
                { if(kv.second != bi) kv.second = MIXED_BUCKET; });
          }
      });

    // Build the finite-cell list once (mesh is read-only during this preprocessing).
    std::vector<typename C3t3::Triangulation::Cell_handle> all_cells;
    {
      const auto& tr_ro = c3t3.triangulation();
      all_cells.reserve(tr_ro.number_of_finite_cells() + 64);
      for(auto cit = tr_ro.finite_cells_begin(); cit != tr_ro.finite_cells_end(); ++cit)
        all_cells.push_back(cit);
    }

    // Pass 2: a cell whose tagged vertices span >=2 buckets (any MIXED, or two
    // different single-bucket tags) is shared; flag all 4 of its vertices
    // boundary-touching. Untagged vertices are ignored (own no element -> no lock).
    boost::concurrent_flat_map<Vh, char, boost::hash<Vh>> boundary_vertex;
    boundary_vertex.reserve(eps.size());
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, all_cells.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        for(std::size_t ci = r.begin(); ci != r.end(); ++ci)
        {
          const auto c = all_cells[ci];
          std::size_t first = 0; bool have = false, shared = false;
          for(int k = 0; k < 4; ++k)
          {
            std::size_t b = 0;
            const bool tagged = vertex_bucket.cvisit(c->vertex(k),
              [&b](const std::pair<const Vh, std::size_t>& kv){ b = kv.second; }) > 0;
            if(!tagged) continue;
            if(b == MIXED_BUCKET || (have && b != first)) { shared = true; break; }
            first = b; have = true;
          }
          if(shared)
            for(int k = 0; k < 4; ++k)
              boundary_vertex.insert_or_assign(c->vertex(k), char(1));
        }
      });

    // Pass 3: process buckets concurrently; each bucket sequential on its worker.
    // INTERIOR elements skip lock_zone (lock-free); BOUNDARY elements lock as before.
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, buckets.size()),
      [&](const tbb::blocked_range<std::size_t>& r)
      {
        boost::container::small_vector<Vh, 2> lv;
        for(std::size_t bi = r.begin(); bi != r.end(); ++bi)
        {
          for(std::size_t i = buckets[bi].first; i < buckets[bi].second; ++i)
          {
            const ElementType& element = eps[i].e;
            lv.clear();
            op.locked_vertices(element, lv);
            bool interior = !lv.empty();
            for(const Vh& v : lv)
              if(boundary_vertex.contains(v)) { interior = false; break; }

            if(interior)
            {
              op.execute_operation(element, c3t3);
#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
              num_elided_locks++;
#endif
            }
            else
            {
              while(!op.lock_zone(element, c3t3))
              {
                c3t3.triangulation().unlock_all_elements();
                std::this_thread::yield();
#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
                num_failed_locks++;
#endif
              }
              op.execute_operation(element, c3t3);
#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
              num_successful_locks++;
#endif
              c3t3.triangulation().unlock_all_elements();
            }
          }
        }
      });

#ifdef CGAL_TETRAHEDRAL_REMESHING_WRITE_LOCK_STATS
    {
      const char* csv_path = "lock_stats.csv";
      bool need_header = true;
      {
        std::ifstream ifs(csv_path);
        need_header = !ifs.good() || (ifs.peek() == std::ifstream::traits_type::eof());
      }
      std::ofstream ofs(csv_path, std::ios::app);
      if(need_header) {
        ofs << "operation,num_successful_locks,num_failed_locks,lock_success_rate(%),num_elided_locks" << std::endl;
      }
      const size_t attempts = static_cast<size_t>(num_successful_locks + num_failed_locks);
      const double lock_success_rate =
          (attempts == 0) ? 0.0 : 100.0 * static_cast<double>(num_successful_locks) / static_cast<double>(attempts);
      ofs << op.operation_name() << "," << num_successful_locks << "," << num_failed_locks << "," << lock_success_rate
          << "," << num_elided_locks << std::endl;
    }
#endif
#endif // concurrent
    return true;
  }

public:
  bool apply_operation_on_elements(std::vector<ElementType>& elements, Operation& op, C3t3& c3t3) override {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Executing operation in parallel: " << op.operation_name() << std::endl;
#endif

    if(elements.empty()) {
      return false;
    }

    // Choose processing strategy based on operation requirements
    if constexpr (Operation::ordered_processing)
      return apply_ordered_processing(elements, op, c3t3, c3t3.bbox());
    else
      return apply_unordered_processing(elements, op, c3t3, c3t3.bbox());
  }
}; // class ElementaryOperationExecution

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H