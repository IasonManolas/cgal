// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_TETRAHEDRAL_REMESHING_IMPL_H
#define CGAL_TETRAHEDRAL_REMESHING_IMPL_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#endif

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/Elementary_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/split_long_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/collapse_short_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/flip_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/smooth_vertices.h>
#include <CGAL/Tetrahedral_remeshing/internal/peel_slivers.h>
#include <CGAL/Tetrahedral_remeshing/internal/property_maps.h>
#include <CGAL/Tetrahedral_remeshing/internal/Parallel_tuning.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/compute_c3t3_statistics.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <optional>
#include <type_traits>
#include <boost/container/small_vector.hpp>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{

class Default_remeshing_visitor
{
public:
  template<typename Tr>
  void before_split(const Tr& /* tr */, const typename Tr::Edge& /* e */) {}
  template<typename Tr>
  void after_split(const Tr& /* tr */, const typename Tr::Vertex_handle /* new_v */) {}

  template<typename CellHandleOld, typename CellHandleNew>
  void after_add_cell(CellHandleOld /* co */, CellHandleNew /* cn */) const {}

  template<typename CellHandle>
  void before_flip(const CellHandle /* c */) {}
  template<typename CellHandle>
  void after_flip(CellHandle /* c */) {}
};



template<typename Triangulation
         , typename SizingFunction
         , typename VertexIsConstrainedMap
         , typename EdgeIsConstrainedMap
         , typename FacetIsConstrainedMap
         , typename CellSelector
         , typename Visitor
         , typename CornerIndex = int
         , typename CurveIndex = int
         >
class Adaptive_remesher
{
  typedef Triangulation Tr;
  typedef typename Tr::Geom_traits::FT FT;

  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex> C3t3;

  typedef typename C3t3::Cell_handle         Cell_handle;
  typedef typename C3t3::Vertex_handle       Vertex_handle;
  typedef typename C3t3::Edge                Edge;
  typedef typename C3t3::Subdomain_index     Subdomain_index;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename C3t3::Curve_index         Curve_index;
  typedef typename C3t3::Corner_index        Corner_index;

  typedef Vertex_smoothing_context<C3t3, SizingFunction, CellSelector> SmoothingContext;

  typedef typename Tr::Concurrency_tag Concurrency_tag;

  // The executor to run an elementary operation with. Sequential unless the
  // triangulation was instantiated with CGAL::Parallel_tag.
  template<typename Operation>
  using Executor = Elementary_operation_executor<Operation, Concurrency_tag>;

private:
  C3t3 m_c3t3;
  const SizingFunction& m_sizing;
  const bool m_protect_boundaries;
  CellSelector m_cell_selector;
  Visitor& m_visitor;
  std::shared_ptr<SmoothingContext> m_smoothing_context;//built from the initial surface

  C3t3* m_c3t3_pbackup;
  std::vector<Vertex_handle> m_far_points;
  Triangulation* m_tr_pbackup; //backup to re-swap triangulations when done

#ifdef CGAL_LINKED_WITH_TBB
  // Owned by the remesher, and handed to the triangulation for the duration of
  // the remeshing, so that the grid matches this triangulation's bounding box
  // and dies with it. Empty under Sequential_tag.
  std::optional<typename Tr::Lock_data_structure> m_lock_ds;
#endif

  /**
  * Cells per axis in the spatial lock grid. 8 is the Mesh_3 default, and is
  * very coarse here: the grid spans the whole bounding box, so a zone covers
  * a large part of the triangulation and two threads working far apart can
  * still contend. A finer grid costs more grid cells to lock per zone and
  * makes a moving vertex cross cells more often, which is why it is a
  * measurable trade rather than an obvious win.
  * `CGAL_TR_LOCK_GRID_SIZE` selects it at run time so both arms live in one
  * binary (POLICY 0.2).
  */
  static int lock_grid_size()
  {
    static const int n = []
      {
        const char* const e = std::getenv("CGAL_TR_LOCK_GRID_SIZE");
        const int v = (e == nullptr) ? 0 : std::atoi(e);
        return (v > 0) ? v : 8;
      }();
    return n;
  }

  /**
  * `CGAL_TR_LOCK_GRID_PER_STAR=k` (A5) sizes the grid from mesh density
  * instead of a constant.
  *
  * The constant is the problem. 8 cells per axis is the Mesh_3 default and is
  * very coarse here; 128 was measured (R3) and lost 2.18%, because a two-ring
  * then spans many grid cells, so both the number of locks per zone and the
  * chance that a moving vertex leaves its cell explode. Neither number is a
  * function of the mesh, and the region between them was never sampled.
  *
  * What should be roughly constant is not cells-per-axis but how many
  * star-sized neighbourhoods fit in a grid cell. Choose N so that a grid cell
  * holds about `k` of them, with a mean star of ~24 cells:
  *
  *     N = cbrt(finite_cells / (k * 24)),  clamped to [4, 64]
  *
  * At 388k cells and k=8 that is N ~ 12; at 1.2M cells, N ~ 17. Both sit in
  * the region between the two constants that have been tried.
  */
  int chosen_lock_grid_size() const
  {
    const int per_star = Tetrahedral_remeshing::internal::Parallel_tuning
                           ::get().lock_grid_per_star;
    if (per_star <= 0)
      return lock_grid_size();

    const double target = static_cast<double>(per_star) * 24.0;
    const double cells = static_cast<double>(
      m_c3t3.triangulation().number_of_finite_cells());
    const int n = static_cast<int>(std::cbrt((std::max)(1.0, cells / target)));
    return std::clamp(n, 4, 64);
  }

  /**
  * Gives the triangulation the lock grid the parallel executors need. Called
  * once the c3t3 is in place, since the grid is sized from its bounding box.
  */
  void init_lock_data_structure()
  {
#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (std::is_convertible_v<Concurrency_tag, CGAL::Parallel_tag>)
    {
      if (m_c3t3.triangulation().get_lock_data_structure() == nullptr)
      {
        m_lock_ds.emplace(m_c3t3.bbox(), chosen_lock_grid_size());
        m_c3t3.triangulation().set_lock_data_structure(std::addressof(*m_lock_ds));
      }
    }
#endif
  }

public:
  Adaptive_remesher(Triangulation& tr
                    , const SizingFunction& sizing
                    , const bool protect_boundaries
                    , VertexIsConstrainedMap vcmap
                    , EdgeIsConstrainedMap ecmap
                    , FacetIsConstrainedMap fcmap
                    , bool smooth_constrained_edges
                    , CellSelector cell_selector
                    , Visitor& visitor
                   )
    : m_c3t3()
    , m_sizing(sizing)
    , m_protect_boundaries(protect_boundaries)
    , m_cell_selector(cell_selector)
    , m_visitor(visitor)
    , m_c3t3_pbackup(NULL)
    , m_tr_pbackup(&tr)
  {
    m_c3t3.triangulation().swap(tr);

    init_c3t3(vcmap, ecmap, fcmap);
    init_lock_data_structure();
    m_smoothing_context = std::make_shared<SmoothingContext>(
      m_c3t3, m_sizing, m_cell_selector, m_protect_boundaries, smooth_constrained_edges);

#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "00-init");
    CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(m_c3t3,
      "00-facets_in_complex_after_init.off");
#endif
  }

  Adaptive_remesher(C3t3& c3t3
                    , const SizingFunction& sizing
                    , const bool protect_boundaries
                    , VertexIsConstrainedMap vcmap
                    , EdgeIsConstrainedMap ecmap
                    , FacetIsConstrainedMap fcmap
                    , bool smooth_constrained_edges
                    , CellSelector cell_selector
                    , Visitor& visitor
                   )
    : m_c3t3()
    , m_sizing(sizing)
    , m_protect_boundaries(protect_boundaries)
    , m_cell_selector(cell_selector)
    , m_visitor(visitor)
    , m_c3t3_pbackup(&c3t3)
    , m_tr_pbackup(NULL)
  {
    m_c3t3.swap(c3t3);

    init_c3t3(vcmap, ecmap, fcmap);
    init_lock_data_structure();
    m_smoothing_context = std::make_shared<SmoothingContext>(
      m_c3t3, m_sizing, m_cell_selector, m_protect_boundaries, smooth_constrained_edges);

#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "00-init");
    CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(m_c3t3,
      "00-facets_in_complex_after_init.off");
#endif
  }

  bool input_is_c3t3() const
  {
    return m_c3t3_pbackup != NULL;
  }

  /**
  * C1 -- one traversal of the finite edges instead of three.
  *
  * Each iteration currently enumerates the edge set three times:
  *
  *     if (!resolution_reached())   // walk 1: stop at the first bad edge
  *     { split();                   // walk 2: collect the too-long edges
  *       collapse(); }              // walk 3: collect the too-short edges
  *
  * and each walk re-runs the same per-edge predicates -- is_in_complex,
  * is_boundary, is_too_long, is_too_short -- over the same edges. Worse,
  * resolution_reached() is not an independent question: it is exactly "no edge
  * is too long or too short", so the first walk looks at the very edges the
  * next two go back to find, and throws them away.
  *
  * MEASURED AND CORRECTED. Fusing all THREE was tried first and is WRONG:
  * split runs between the second and third traversal and changes the mesh, so
  * a collapse candidate list computed before split names edges split has since
  * destroyed and misses the ones it created. It scored -8.55% wall and, more
  * to the point, BLOCKED on correctness -- 1.92% more cells than the baseline,
  * with three Tier-A configs past the per-config quality limit. The three
  * traversals are not redundant with each other; only the first two are, and
  * that is the whole extent of what can be fused.
  *
  * So this fuses resolution_reached() and split's collection, which do see the
  * same mesh, into ONE parallel cell scan. Collapse keeps its own scan. Three
  * traversals become two.
  *
  * The resolution test is computed independently of can_be_split/
  * can_be_collapsed rather than inferred from whether the two lists are empty:
  * those predicates are strictly narrower (an edge can be too long and not
  * splittable), so inferring it would change the loop's termination.
  *
  * This also makes R9 -- the parallel resolution scan, rejected at -0.49%
  * because a parallel scan cannot early-exit as promptly as the serial one --
  * free rather than negative: the scan is no longer optional, so the early
  * exit it gave up was never worth anything here.
  */
  using Edge_with_length = std::pair<Edge, FT>;
  struct Fused_edges
  {
    std::vector<Edge_with_length> too_long;
    bool resolution_reached = true;
  };
  Fused_edges m_fused;
  bool m_fused_valid = false;

  void run_fused_edge_pass()
  {
    m_fused.too_long.clear();
    m_fused.resolution_reached = true;
    m_fused_valid = false;

#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (std::is_convertible_v<Concurrency_tag, CGAL::Parallel_tag>)
    {
      struct Rec { Edge e; FT sqlen; unsigned char kind; }; // 1 long, 2 short
      std::atomic<bool> resolved{true};
      const Tr& tr = this->tr();

      const std::vector<Rec> found = parallel_collect_finite_edges<Rec>(
        tr,
        [&](const Edge& e, std::vector<Rec>& out)
        {
          const bool boundary =
            m_c3t3.is_in_complex(e) || is_boundary(m_c3t3, e, m_cell_selector);
          if (!(m_protect_boundaries && boundary))
          {
            if (is_too_long(e, boundary, m_sizing, m_c3t3, m_cell_selector)
             || is_too_short(e, boundary, m_sizing, m_c3t3, m_cell_selector))
              resolved.store(false, std::memory_order_relaxed);
          }

          {
            auto [splittable, b] = can_be_split(e, m_c3t3, m_protect_boundaries, m_cell_selector);
            if (splittable)
            {
              const std::optional<FT> sqlen
                = is_too_long(e, b, m_sizing, m_c3t3, m_cell_selector);
              if (sqlen != std::nullopt)
                out.push_back(Rec{e, sqlen.value(), 1});
            }
          }
        });

      for (const Rec& r : found)
        m_fused.too_long.emplace_back(r.e, r.sqlen);
      m_fused.resolution_reached = resolved.load();
      m_fused_valid = true;
    }
#endif
  }

  bool fused_enabled() const
  {
    if constexpr (std::is_convertible_v<Concurrency_tag, CGAL::Parallel_tag>)
      return Tetrahedral_remeshing::internal::Parallel_tuning::get().fused_edge_pass;
    else
      return false;
  }

  void split()
  {
    CGAL_assertion(check_vertex_dimensions());
    typedef Edge_split_operation<C3t3, SizingFunction, CellSelector, Visitor> EdgeSplitOp;
    EdgeSplitOp split_op(m_sizing, m_cell_selector, m_protect_boundaries, m_visitor);
    if (m_fused_valid)
      split_op.set_precollected(&m_fused.too_long);
    Executor<EdgeSplitOp> executor;
    executor.execute(split_op, m_c3t3);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
    CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(m_c3t3,
      "1-facets_in_complex_after_split.off");
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "1-c3t3_vertices_after_split");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
    CGAL::Tetrahedral_remeshing::debug::dump_edges_in_complex(m_c3t3,
      "1-edges_in_complex_after_split.polylines.txt");
    const double mdh = CGAL::Tetrahedral_remeshing::min_dihedral_angle(m_c3t3);
    std::cout << "Min dihedral angle = " << mdh << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "1-split");
#endif
  }

  void collapse()
  {
    CGAL_assertion(check_vertex_dimensions());
    typedef Edge_collapse_operation<C3t3, SizingFunction, CellSelector, Visitor> EdgeCollapseOp;
    EdgeCollapseOp collapse_op(m_sizing, m_cell_selector, m_protect_boundaries, m_visitor);
    Executor<EdgeCollapseOp> executor;
    executor.execute(collapse_op, m_c3t3);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "2-c3t3_vertices_after_collapse");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
    CGAL::Tetrahedral_remeshing::debug::dump_edges_in_complex(m_c3t3,
      "2-edges_in_complex_after_collapse.polylines.txt");
    const double mdh = CGAL::Tetrahedral_remeshing::min_dihedral_angle(m_c3t3);
    std::cout << "Min dihedral angle = " << mdh << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "2-collapse");
#endif
  }

  void flip()
  {
    typedef Internal_edge_flip_operation<C3t3, CellSelector, Visitor> InternalFlipOp;
    typedef Boundary_edge_flip_operation<C3t3, CellSelector, Visitor> BoundaryFlipOp;

    // one incident-cells cache shared by both passes, as in the former flip_edges()
    typename InternalFlipOp::Incident_cells_map inc_cells;

    InternalFlipOp internal_flip_op(m_cell_selector, m_visitor, inc_cells);
    Executor<InternalFlipOp> internal_executor;
    internal_executor.execute(internal_flip_op, m_c3t3);

    if (!m_protect_boundaries)
    {
      BoundaryFlipOp boundary_flip_op(m_cell_selector, m_visitor, inc_cells);
      Executor<BoundaryFlipOp> boundary_executor;
      boundary_executor.execute(boundary_flip_op, m_c3t3);
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "3-c3t3_vertices_after_flip");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
    CGAL::Tetrahedral_remeshing::debug::dump_edges_in_complex(m_c3t3,
      "3-edges_in_complex_after_flip.polylines.txt");
    const double mdh = CGAL::Tetrahedral_remeshing::min_dihedral_angle(m_c3t3);
    std::cout << "Min dihedral angle = " << mdh << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "3-flip");
#endif
  }

  void smooth()
  {
    m_smoothing_context->refresh(m_c3t3);

    // Order matches the former Tetrahedral_remeshing_smoother::smooth_vertices():
    // complex (1D) edges, then surface (2D) vertices, then internal (3D) vertices.
    if (!m_protect_boundaries)
    {
      if (m_smoothing_context->m_smooth_constrained_edges)
      {
        Complex_edge_vertex_smooth_operation<C3t3, SizingFunction, CellSelector> op(m_smoothing_context);
        Executor<decltype(op)> executor;
        executor.execute(op, m_c3t3);
      }
      {
        Surface_vertex_smooth_operation<C3t3, SizingFunction, CellSelector> op(m_smoothing_context);
        Executor<decltype(op)> executor;
        executor.execute(op, m_c3t3);
      }
    }
    {
      Internal_vertex_smooth_operation<C3t3, SizingFunction, CellSelector> op(m_smoothing_context);
      Executor<decltype(op)> executor;
      executor.execute(op, m_c3t3);
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "4-c3t3_vertices_after_smooth");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
    CGAL::Tetrahedral_remeshing::debug::dump_edges_in_complex(m_c3t3,
      "4-edges_in_complex_after_smoothing.polylines.txt");
    const double mdh = CGAL::Tetrahedral_remeshing::min_dihedral_angle(m_c3t3);
    std::cout << "Min dihedral angle = " << mdh << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "4-smooth");
#endif
  }

  /**
  * `CGAL_TR_PARALLEL_RESOLUTION=1` answers the same question with a parallel
  * cell scan and an atomic early exit. Off by default.
  *
  * The serial loop stops at the FIRST edge that is still too long or too
  * short, so when the answer is "not yet" it usually costs very little -- and
  * that is the common case until the last iteration. The parallel arm cannot
  * stop as promptly: threads already inside a block run on. So this trades a
  * cheap early exit for parallelism that only pays on the one iteration where
  * the answer is "yes" and every edge must be examined.
  */
  static bool parallel_resolution_enabled()
  {
    static const bool enabled = []
      {
        const char* const e = std::getenv("CGAL_TR_PARALLEL_RESOLUTION");
        return (e != nullptr) && (std::atoi(e) != 0);
      }();
    return enabled;
  }

  bool resolution_reached()
  {
#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (std::is_convertible_v<Concurrency_tag, CGAL::Parallel_tag>)
    {
      if (parallel_resolution_enabled())
      {
        std::atomic<bool> done{true};
        const std::vector<Edge> edges
          = parallel_collect_finite_edges<Edge>(
              tr(),
              [](const Edge& e, std::vector<Edge>& out) { out.push_back(e); });

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, edges.size()),
          [&](const tbb::blocked_range<std::size_t>& r)
          {
            if (!done.load(std::memory_order_relaxed))
              return; // another thread already found one; do not start a block
            for (std::size_t i = r.begin(); i != r.end(); ++i)
            {
              const Edge& e = edges[i];
              const bool boundary =
                m_c3t3.is_in_complex(e) || is_boundary(m_c3t3, e, m_cell_selector);
              if (m_protect_boundaries && boundary)
                continue;
              if (is_too_long(e, boundary, m_sizing, m_c3t3, m_cell_selector)
               || is_too_short(e, boundary, m_sizing, m_c3t3, m_cell_selector))
              {
                done.store(false, std::memory_order_relaxed);
                return;
              }
            }
          });
        return done.load();
      }
    }
#endif

    for (const Edge& e : tr().finite_edges())
    {
      // skip protected edges
      const bool boundary =
        m_c3t3.is_in_complex(e) || is_boundary(m_c3t3, e, m_cell_selector);
      if (m_protect_boundaries && boundary)
        continue;

      if(  is_too_long(e, boundary, m_sizing, m_c3t3, m_cell_selector)
        || is_too_short(e, boundary, m_sizing, m_c3t3, m_cell_selector))
        return false;
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Resolution reached" << std::endl;
#endif
    return true;
  }

  //peel off slivers
  std::size_t postprocess(const double sliver_angle = 2.)
  {
    if (m_protect_boundaries)
      return 0;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Postprocess...";
    std::cout.flush();
#endif

    const std::size_t nb_peeled
      = CGAL::Tetrahedral_remeshing::peel_slivers(m_c3t3, sliver_angle, m_cell_selector);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "99-postprocess");
#endif

    return nb_peeled;
  }

  void finalize()
  {
#ifdef CGAL_LINKED_WITH_TBB
    // The grid dies with the remesher. Take it back before the triangulation
    // is swapped out, so the caller's does not leave holding a stale pointer.
    if (m_lock_ds.has_value())
      m_c3t3.triangulation().set_lock_data_structure(nullptr);
#endif

    if (input_is_c3t3())
    {
      //reset far points dimension
      for (Vertex_handle v : m_far_points)
        v->set_dimension(-1);
      m_c3t3_pbackup->swap(m_c3t3);
    }
    else
      m_tr_pbackup->swap(m_c3t3.triangulation());
  }

private:
  void init_c3t3(const VertexIsConstrainedMap& vcmap,
                 const EdgeIsConstrainedMap& ecmap,
                 const FacetIsConstrainedMap& fcmap)
  {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    debug_c3t3();
    std::size_t nbc = 0;
    std::size_t nbf = 0;
    std::size_t nbe = 0;
    std::size_t nbv = 0;
#endif
    //update number_of_cells and number_of_facets in c3t3
    m_c3t3.rescan_after_load_of_triangulation();

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "00-c3t3_vertices_before_init_");
#endif

    if (input_is_c3t3())
      backup_far_points();

    const Subdomain_index default_subdomain = default_subdomain_index();

    //tag cells
    for (Cell_handle cit : tr().finite_cell_handles())
    {
      if (get(m_cell_selector, cit))
      {
        const Subdomain_index index = cit->subdomain_index();
        if (m_c3t3.is_in_complex(cit))
          m_c3t3.remove_from_complex(cit);

        const Subdomain_index new_index = (Subdomain_index() != index)
          ? index
          : default_subdomain;
        m_c3t3.add_to_complex(cit, new_index);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ++nbc;
#endif
      }

      for (Vertex_handle vi : tr().vertices(cit))
        set_dimension(vi, 3);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      if (input_is_c3t3() && m_c3t3.is_in_complex(cit))
        ++nbc;
#endif
    }

    //tag facets
    typedef typename Tr::Facet Facet;
    for (const Facet& f : tr().finite_facets())
    {
      const Facet mf = tr().mirror_facet(f);
      const Subdomain_index s1 = f.first->subdomain_index();
      const Subdomain_index s2 = mf.first->subdomain_index();
      if (s1 != s2
          || get(fcmap, f)
          || get(fcmap, mf)
          || (input_is_c3t3() && m_c3t3.is_in_complex(f))
          || (!input_is_c3t3() && f.first->is_facet_on_surface(f.second)))
      {
        Surface_patch_index patch = f.first->surface_patch_index(f.second);
        if(patch == Surface_patch_index())
          make_surface_patch_index(s1, s2, patch);

        if(m_c3t3.is_in_complex(f))
          m_c3t3.remove_from_complex(f);
        m_c3t3.add_to_complex(f, patch);

        for (Vertex_handle vij : tr().vertices(f))
          set_dimension(vij, 2);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ++nbf;
#endif
      }
    }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(m_c3t3, "facets_in_complex.off");
#endif

    //tag edges
    const Curve_index default_curve_id = default_curve_index();
    for (const Edge& e : tr().finite_edges())
    {
      if (get(ecmap, CGAL::Tetrahedral_remeshing::make_vertex_pair(e))
          || get(ecmap, CGAL::Tetrahedral_remeshing::make_inv_vertex_pair(e))
          || (input_is_c3t3() && m_c3t3.is_in_complex(e))
          || nb_incident_subdomains(e, m_c3t3) > 2
          || nb_incident_surface_patches(e, m_c3t3) > 1
          || nb_incident_complex_facets(e, m_c3t3) > 2)//non-manifold edges
      {
        const bool in_complex = m_c3t3.is_in_complex(e);
        typename C3t3::Curve_index curve_id = in_complex
          ? m_c3t3.curve_index(e)
          : default_curve_id;

        if (in_complex)
          m_c3t3.remove_from_complex(e);
        m_c3t3.add_to_complex(e, curve_id);

        for (Vertex_handle v : tr().vertices(e))
          set_dimension(v, 1);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ++nbe;
#endif
      }
    }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL::Tetrahedral_remeshing::debug::dump_edges_in_complex(m_c3t3, "edges_in_complex.polylines.txt");
#endif

    //tag vertices
    Corner_index corner_id = 0;
    for (Vertex_handle vit : tr().finite_vertex_handles())
    {
      boost::container::small_vector<Edge, 3> incident_edges;
      incident_complex_edges(vit, m_c3t3, std::back_inserter(incident_edges));
      if ( incident_edges.size() == 1 //tip or endpoint
        || incident_edges.size() > 2  //corner
        || get(vcmap, vit) // constrained vertex
        || (input_is_c3t3() && m_c3t3.is_in_complex(vit))
        || (incident_edges.size() == 2
            && edges_form_a_sharp_angle(incident_edges, 60, m_c3t3)))
      {
        if (!m_c3t3.is_in_complex(vit))
          m_c3t3.add_to_complex(vit, ++corner_id);

        set_dimension(vit, 0);
        vit->set_index(corner_id);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ++nbv;
#endif
      }
    }

    // set all indices depending on underlying dimension
    for (Vertex_handle v : tr().finite_vertex_handles())
      set_index(v, m_c3t3);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    std::cout << "C3t3 ready :" << std::endl;
    std::cout << "\t cells    = " << nbc << std::endl;
    std::cout << "\t facets   = " << nbf << std::endl;
    std::cout << "\t edges    = " << nbe << std::endl;
    std::cout << "\t vertices = " << nbv << std::endl;
    const double mdh = CGAL::Tetrahedral_remeshing::min_dihedral_angle(m_c3t3);
    std::cout << "\t Min dihedral angle = " << mdh << std::endl;

    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "0-c3t3_vertices_after_init_");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
    CGAL::Tetrahedral_remeshing::debug::count_far_points(m_c3t3);
    CGAL::Tetrahedral_remeshing::debug::dump_edges_in_complex(m_c3t3,
      "0-edges_in_complex_after_init.polylines.txt");
#endif
  }

  void backup_far_points()
  {
    for (Vertex_handle v : tr().finite_vertex_handles())
    {
      if (v->in_dimension() == -1)
        m_far_points.push_back(v);
    }
  }

private:
  bool dimension_is_modifiable(const Vertex_handle& v, const int new_dim) const
  {
    const int vdim = v->in_dimension();
    // feature edges and tip/endpoints vertices are kept
    switch (vdim)
    {
    case -1: return (new_dim == 3);//far points
    case 3 : return true;
    case 2 : return true;//surface vertices may not be part of a triangulation surface
                        // in this case, we want to be able to set it
    case 1 : return new_dim == 0; //features can be modified to corners
    case 0 : return false;// corners remain corners
    default: break;
    }
    CGAL_assertion(false);
    return true;
  }

  void set_dimension(Vertex_handle v, const int new_dim)
  {
    if (dimension_is_modifiable(v, new_dim))
      v->set_dimension(new_dim);
  }

  bool check_vertex_dimensions()
  {
    for (Vertex_handle vit : tr().finite_vertex_handles())
    {
      // dimension is -1 for Mesh_3 "far points"
      // for other vertices, it is in [0; 3]
      if (vit->in_dimension() < -1 || vit->in_dimension() > 3)
        return false;
    }
    return true;
  }
  void debug_c3t3()
  {
    CGAL_assertion_code(for (typename Tr::Facet f : tr().finite_facets()))
    {
      CGAL_assertion_code(typename Tr::Facet mf = tr().mirror_facet(f));
      CGAL_assertion(m_c3t3.is_in_complex(f) == m_c3t3.is_in_complex(mf));
    }
  }

  template<typename PatchIndex>
  void make_surface_patch_index(const Subdomain_index& s1,
                                const Subdomain_index& s2,
                                PatchIndex& patch)
  {
    patch = (s1 < s2) ? (s1 * 1000 + s2) : (s2 * 1000 + s1);
  }

  void make_surface_patch_index(const Subdomain_index& s1,
                                const Subdomain_index& s2,
                                std::pair<Subdomain_index, Subdomain_index>& patch)
  {
    patch = (s1 < s2) ? std::make_pair(s1, s2) : std::make_pair(s2, s1);
  }

  Subdomain_index max_subdomain_index() const
  {
    Subdomain_index max_index
      = (std::numeric_limits<Subdomain_index>::min)();
    for (Cell_handle cit : tr().finite_cell_handles())
    {
      const Subdomain_index cid = cit->subdomain_index();
      if (cid > max_index && cid != Subdomain_index())
        max_index = cid;
    }
    return max_index;
  }

  Subdomain_index default_subdomain_index() const
  {
    return max_subdomain_index() + 1;
  }

  Curve_index max_curve_index() const
  {
    Curve_index max_index = (std::numeric_limits<Curve_index>::min)();
    for (const Edge& e : m_c3t3.edges_in_complex())
    {
      const Curve_index cid = m_c3t3.curve_index(e);
      if (cid > max_index)
        max_index = cid;
    }
    if (max_index == (std::numeric_limits<Curve_index>::min)())
      return 0;
    else
      return max_index + 1;
  }

  Curve_index default_curve_index() const
  {
    return max_curve_index() + 1;
  }

public:
  Tr& tr()
  {
    return m_c3t3.triangulation();
  }
  const Tr& tr() const
  {
    return m_c3t3.triangulation();
  }

  void remesh(const std::size_t& max_it,
              const std::size_t& nb_extra_iterations)
  {
    std::size_t it_nb = 0;
    while (it_nb < max_it)
    {
      ++it_nb;
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " #" << std::endl;
#endif
      // D2: the shared spatial subdivision is rebuilt once per iteration, so
      // every unordered operation in the iteration uses the same regions.
#ifdef CGAL_LINKED_WITH_TBB
      if constexpr (std::is_convertible_v<Concurrency_tag, CGAL::Parallel_tag>)
        if (Tetrahedral_remeshing::internal::Parallel_tuning::get().reuse_partition)
          Tetrahedral_remeshing::internal::Shared_kd_partition::get().reset();
#endif
      bool resolved;
      if (fused_enabled())
      {
        run_fused_edge_pass();
        resolved = m_fused.resolution_reached;
      }
      else
      {
        m_fused_valid = false;
        resolved = resolution_reached();
      }
      if (!resolved)
      {
        split();
        collapse();
      }
      m_fused_valid = false;
      flip();
      smooth();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " done : "
                << tr().number_of_vertices()
                << " vertices #" << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
      std::ostringstream ossi;
      ossi << "statistics_" << it_nb << ".txt";
      Tetrahedral_remeshing::internal::compute_statistics(
        tr(), m_cell_selector, ossi.str().c_str());
      std::ostringstream oss_it;
      oss_it << "iteration_" << it_nb;
      Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, oss_it.str().c_str());
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
#endif
    }

    m_smoothing_context->start_flip_smooth_steps(m_c3t3);
    while (it_nb < max_it + nb_extra_iterations)
    {
      ++it_nb;

#ifdef CGAL_LINKED_WITH_TBB
      if constexpr (std::is_convertible_v<Concurrency_tag, CGAL::Parallel_tag>)
        if (Tetrahedral_remeshing::internal::Parallel_tuning::get().reuse_partition)
          Tetrahedral_remeshing::internal::Shared_kd_partition::get().reset();
#endif
      flip();
      smooth();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " (flip and smooth only) done : "
                << tr().number_of_vertices()
                << " vertices #" << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
      std::ostringstream ossi;
      ossi << "statistics_" << it_nb << ".txt";
      Tetrahedral_remeshing::internal::compute_statistics(
        tr(),  m_cell_selector, ossi.str().c_str());
#endif
    }

    postprocess(); //peel off boundary slivers

    finalize();
    //Warning : triangulation() is now empty
  }

};//end class Adaptive_remesher


template<typename Triangulation,
         typename SizingFunction,
         typename NamedParameters,
         typename CornerIndex = int,
         typename CurveIndex = int>
struct Adaptive_remesher_type_generator
{
  using Tr = Triangulation;

  using Default_Selection_functor = All_cells_selected<Tr>;
  using SelectionFunctor = typename internal_np::Lookup_named_param_def<
    internal_np::cell_selector_t,
    NamedParameters,
    Default_Selection_functor//default
  >::type;

  using Vertex_handle = typename Tr::Vertex_handle;
  using Default_VCMap = Static_boolean_property_map<Vertex_handle, false>;
  using VCMap = typename internal_np::Lookup_named_param_def<
    internal_np::vertex_is_constrained_t,
    NamedParameters,
    Default_VCMap//default
  >::type;

  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using Default_ECMap = Static_boolean_property_map<Edge_vv, false>;
  using ECMap = typename internal_np::Lookup_named_param_def<
    internal_np::edge_is_constrained_t,
    NamedParameters,
    Default_ECMap//default
  >::type;

  using Facet = typename Tr::Facet;
  using Default_FCMap = Static_boolean_property_map<Facet, false>;
  using FCMap = typename internal_np::Lookup_named_param_def<
    internal_np::facet_is_constrained_t,
    NamedParameters,
    Default_FCMap//default
  >::type;

  using Default_Visitor = Default_remeshing_visitor;
  using Visitor = typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters,
    Default_Visitor//default
  >::type;

  using type = Adaptive_remesher<
    Tr, SizingFunction, VCMap, ECMap, FCMap, SelectionFunctor, Visitor>;
};

}//end namespace internal
}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif //CGAL_TETRAHEDRAL_REMESHING_IMPL_H
