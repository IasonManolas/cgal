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

#ifndef CGAL_INTERNAL_SPLIT_LONG_EDGES_H
#define CGAL_INTERNAL_SPLIT_LONG_EDGES_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_sort.h>
#endif
#include <CGAL/Tetrahedral_remeshing/internal/Parallel_tuning.h>

#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include <CGAL/Tetrahedral_remeshing/internal/Elementary_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/MVLZ_probe.h>

#include <unordered_map>
#include <functional>
#include <utility>
#include <optional>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{

template<typename C3t3>
bool positive_orientation_after_edge_split(const typename C3t3::Edge& e,
                                           const typename C3t3::Cell_handle circ,
                                           const typename C3t3::Triangulation::Geom_traits::Point_3& steiner,
                                           const C3t3&)
{
  using Point = typename C3t3::Triangulation::Geom_traits::Point_3;

  const auto v1 = e.first->vertex(e.second);
  const auto v2 = e.first->vertex(e.third);

  std::array<Point, 4> pts = {point(circ->vertex(0)->point()),
                              point(circ->vertex(1)->point()),
                              point(circ->vertex(2)->point()),
                              point(circ->vertex(3)->point())};
  // 1st half-cell
  const int i1 = circ->index(v1);
  const Point p1 = pts[i1];
  pts[i1] = steiner;
  if(CGAL::orientation(pts[0], pts[1], pts[2], pts[3]) != CGAL::POSITIVE)
    return false;

  // 2nd half-cell
  pts[i1] = p1;
  pts[circ->index(v2)] = steiner;
  if(CGAL::orientation(pts[0], pts[1], pts[2], pts[3]) != CGAL::POSITIVE)
    return false;

  return true;
}

template <typename C3t3>
std::optional<typename C3t3::Triangulation::Geom_traits::Point_3>
construct_steiner_point(const typename C3t3::Edge& e,
                        const C3t3& c3t3)
{
  using Cell_circulator = typename C3t3::Triangulation::Cell_circulator;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;
  using Point = typename C3t3::Triangulation::Geom_traits::Point_3;
  using FT = typename C3t3::Triangulation::Geom_traits::FT;

  const auto& gt = c3t3.triangulation().geom_traits();
  const auto& tr = c3t3.triangulation();
  const auto& p1 = point(e.first->vertex(e.second)->point());
  const auto& p2 = point(e.first->vertex(e.third)->point());
  const auto vec = gt.construct_vector_3_object()(p1, p2);

  const std::array<FT, 6> coeff = {0.33, 0.66,    //1/3 and 2/3
                                   0.3, 0.7,      // 0.5 +/- 0.2
                                   0.25, 0.75};   // 0.5 +/- 0.25

  std::size_t attempt_id = 0;
  while(attempt_id < coeff.size())
  {
    Point steiner = gt.construct_translated_point_3_object()(
        p1, gt.construct_scaled_vector_3_object()(vec, coeff[attempt_id]));
    ++attempt_id;

    bool steiner_successful = true;
    Cell_circulator circ = tr.incident_cells(e);
    Cell_circulator end = circ;
    do
    {
      Cell_handle c = circ;
      if(!positive_orientation_after_edge_split(e, c, steiner, c3t3))
      {
        steiner_successful = false;
        break;
      }
    } while(++circ != end);

    if(steiner_successful)
      return steiner;
  }

  return std::nullopt;
}

template<typename C3t3, typename CellSelector>
typename C3t3::Vertex_handle split_edge(const typename C3t3::Edge& e,
                                        CellSelector cell_selector,
                                        C3t3& c3t3)
{
  typedef typename C3t3::Triangulation       Tr;
  typedef typename C3t3::Subdomain_index     Subdomain_index;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename C3t3::Curve_index         Curve_index;
  typedef typename Tr::Geom_traits::Point_3 Point;
  typedef typename Tr::Facet                Facet;
  typedef typename Tr::Vertex_handle        Vertex_handle;
  typedef typename Tr::Cell_handle          Cell_handle;
  typedef typename Tr::Cell_circulator      Cell_circulator;

  Tr& tr = c3t3.triangulation();
  const Vertex_handle v1 = e.first->vertex(e.second);
  const Vertex_handle v2 = e.first->vertex(e.third);

  Point m = tr.geom_traits().construct_midpoint_3_object()
    (point(v1->point()), point(v2->point()));

  //backup subdomain info of incident cells before making changes
  short dimension = 0;
  if(c3t3.is_in_complex(e))
    dimension = 1;
  else
  {
    const std::size_t nb_patches = nb_incident_surface_patches(e, c3t3);
    if(nb_patches == 1)
      dimension = 2;
    else if(nb_patches == 0)
      dimension = 3;
    else
      CGAL_assertion(false);//e should be in complex
  }
  CGAL_assertion(dimension > 0);

  // remove complex edge before splitting
  const Curve_index curve_index = (dimension == 1) ? c3t3.curve_index(e) : Curve_index();

  struct Cell_info {
    Subdomain_index subdomain_index_;
    bool selected_;
  };
  struct Facet_info {
    Vertex_handle opp_vertex_;
    Surface_patch_index patch_index_;
  };
  boost::unordered_map<Facet, Cell_info, boost::hash<Facet>> cells_info;
  boost::unordered_map<Facet, Facet_info, boost::hash<Facet>> facets_info;

  // check orientation and collect incident cells to avoid circulating twice
  bool steiner_point_found = false;
  boost::container::small_vector<Cell_handle, 30> inc_cells;
  Cell_circulator circ = tr.incident_cells(e);
  Cell_circulator end = circ;
  do
  {
    inc_cells.push_back(circ);
    if (tr.is_infinite(circ) || steiner_point_found)
    {
      ++circ;
      continue;
    }

    const Cell_handle c = circ;
    if(!positive_orientation_after_edge_split(e, c, m, c3t3))
    {
      const std::optional<Point> steiner = construct_steiner_point(e, c3t3);
      if (steiner != std::nullopt)
      {
        m = *steiner;
        steiner_point_found = true;
      }
      else
        return Vertex_handle();
    }
    ++circ;
  }
  while (circ != end);

  if (dimension == 1)
    c3t3.remove_from_complex(e);

  for(Cell_handle c : inc_cells)
  {
    const int index_v1 = c->index(v1);
    const int index_v2 = c->index(v2);

    //keys are the opposite facets to the ones not containing e,
    //because they will not be modified
    const Subdomain_index subdomain = c3t3.subdomain_index(c);
    const bool selected = get(cell_selector, c);
    const Facet opp_facet1 = tr.mirror_facet(Facet(c, index_v1));
    const Facet opp_facet2 = tr.mirror_facet(Facet(c, index_v2));

    // volume data
    cells_info.insert(std::make_pair(opp_facet1, Cell_info{subdomain, selected}));
    cells_info.insert(std::make_pair(opp_facet2, Cell_info{subdomain, selected}));
    treat_before_delete(c, cell_selector, c3t3);

    // surface data for facets of the cells to be split
    const int findex = CGAL::Triangulation_utils_3::next_around_edge(index_v1, index_v2);
    Surface_patch_index patch = c3t3.surface_patch_index(c, findex);
    Vertex_handle opp_vertex = c->vertex(findex);
    facets_info.insert(std::make_pair(opp_facet1, Facet_info{opp_vertex, patch}));
    facets_info.insert(std::make_pair(opp_facet2, Facet_info{opp_vertex, patch}));

    if(c3t3.is_in_complex(c, findex))
      c3t3.remove_from_complex(c, findex);
  }

  // insert midpoint
  CGAL_TR_PROBE_CELL_WRITE(tr, e.first, "split: insert_in_edge");
  Vertex_handle new_v = tr.tds().insert_in_edge(e);
  new_v->set_point(typename Tr::Point(m));
  new_v->set_dimension(dimension);

  // update c3t3 with subdomain and surface patch indices
  std::vector<Cell_handle> new_cells;
  tr.incident_cells(new_v, std::back_inserter(new_cells));
  for (Cell_handle new_cell : new_cells)
  {
    const Facet fi(new_cell, new_cell->index(new_v));
    const Facet mfi = tr.mirror_facet(fi);

    //get subdomain info back
    CGAL_assertion(cells_info.find(mfi) != cells_info.end());
    Cell_info c_info = cells_info.at(mfi);
    treat_new_cell(new_cell, c_info.subdomain_index_,
                   cell_selector, c_info.selected_, c3t3);

    // get surface info back
    CGAL_assertion(facets_info.find(mfi) != facets_info.end());
    const Facet_info v_and_opp_patch = facets_info.at(mfi);

    // facet opposite to new_v (status wrt c3t3 is unchanged)
    new_cell->set_surface_patch_index(new_cell->index(new_v),
                                      mfi.first->surface_patch_index(mfi.second));

    // new half-facet (added or not to c3t3 depending on the stored surface patch index)
    if (Surface_patch_index() == v_and_opp_patch.patch_index_)
      new_cell->set_surface_patch_index(new_cell->index(v_and_opp_patch.opp_vertex_),
                                        Surface_patch_index());
    else
      c3t3.add_to_complex(new_cell,
                          new_cell->index(v_and_opp_patch.opp_vertex_),
                          v_and_opp_patch.patch_index_);

    // newly created internal facet
    for (int i = 0; i < 4; ++i)
    {
      const Vertex_handle vi = new_cell->vertex(i);
      if (vi == v1 || vi == v2)
      {
        new_cell->set_surface_patch_index(i, Surface_patch_index());
        break;
      }
    }

    //the 4th facet (new_v, v_and_opp_patch.first, v1 or v2)
    // will have its patch tagged from the other side, if needed
  }

  // re-insert complex sub-edges
  if (dimension == 1)
  {
    c3t3.add_to_complex(new_v, v1, curve_index);
    c3t3.add_to_complex(new_v, v2, curve_index);
  }

  set_index(new_v, c3t3);

  return new_v;
}

/**
* returns [can_be_split, is_on_boundary]
*/
template<typename C3T3, typename CellSelector>
auto can_be_split(const typename C3T3::Edge& e,
                  const C3T3& c3t3,
                  const bool protect_boundaries,
                  const CellSelector& cell_selector)
{
  struct Splittable
  {
    bool can_be_split;
    bool on_boundary;
  };

  if (is_outside(e, c3t3, cell_selector))
    return Splittable{false, false};

  const bool boundary = c3t3.is_in_complex(e)
                     || is_boundary(c3t3, e, cell_selector);

  if (protect_boundaries)
  {
    if (boundary)
      return Splittable{false, boundary};

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    if (!is_internal(e, c3t3, cell_selector))
    {
      std::cerr << "e is not inside!?" << std::endl;
      typename C3T3::Vertex_handle v1 = e.first->vertex(e.second);
      typename C3T3::Vertex_handle v2 = e.first->vertex(e.third);
      std::cerr << v1->point() << " " << v2->point() << std::endl;
    }
#endif

    CGAL_assertion(is_internal(e, c3t3, cell_selector));
    return Splittable{true, boundary};
  }
  else
  {
    return Splittable{is_selected(e, c3t3.triangulation(), cell_selector), boundary};
  }
}



template<typename C3t3,
         typename SizingFunction,
         typename CellSelector,
         typename Visitor>
class Edge_split_operation
    : public Elementary_operation<C3t3,
                                 std::pair<typename C3t3::Triangulation::Vertex_handle,
                                           typename C3t3::Triangulation::Vertex_handle>,
                                 std::vector<std::pair<typename C3t3::Triangulation::Vertex_handle,
                                                       typename C3t3::Triangulation::Vertex_handle>>>
{
public:
  using Tr = typename C3t3::Triangulation;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Cell_handle = typename Tr::Cell_handle;
  using Edge = typename Tr::Edge;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using FT = typename Tr::Geom_traits::FT;

  // Candidates are stored as vertex pairs, captured at collection time. A raw
  // Edge (Cell_handle, i, j) would go stale: each split destroys and recycles
  // cells, so by the time the executor reaches a later candidate its Cell_handle
  // may point at a different cell. Vertices are never removed by a split, so the
  // vertex pair stays valid and is re-resolved to the current edge via is_edge().
  using Long_edges = std::vector<Edge_vv>;
  using Base_operation = Elementary_operation<C3t3, Edge_vv, Long_edges>;
  using Element_type = typename Base_operation::Element_type;
  static_assert(std::is_same_v<Element_type, Edge_vv>, "Element_type must be Edge_vv");
  using ElementSource = typename Base_operation::Element_range;

private:
  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  Visitor& m_visitor;

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
  mutable std::ofstream m_can_be_split_ofs;
  mutable std::ofstream m_split_failed_ofs;
  mutable std::ofstream m_midpoints_ofs;
#endif

public:
  Edge_split_operation(const SizingFunction& sizing,
                     const CellSelector& cell_selector,
                     const bool protect_boundaries,
                     Visitor& visitor)
      : m_sizing(sizing)
      , m_cell_selector(cell_selector)
      , m_protect_boundaries(protect_boundaries)
      , m_visitor(visitor) {}

  using Edge_with_length = std::pair<Edge, FT>;

  /** C1: candidates already collected by the fused edge pass, or nullptr. */
  void set_precollected(const std::vector<Edge_with_length>* p) { m_precollected = p; }

  ElementSource get_elements(const C3t3& c3t3) const override
  {
    struct Long_edge_with_length
    {
      Edge edge;
      FT sqlength;
    };
    std::vector<Long_edge_with_length> long_edges_with_lengths;
    const Tr& tr = c3t3.triangulation();

    // The test applied to each edge, shared by the serial walk and the
    // parallel cell scan so the two collect exactly the same set.
    const auto keep = [&](const Edge& e, std::vector<Long_edge_with_length>& out)
    {
      auto [splittable, boundary] = can_be_split(e, c3t3, m_protect_boundaries, m_cell_selector);
      if (!splittable)
        return;

      const std::optional<FT> sqlen = is_too_long(e, boundary, m_sizing, c3t3, m_cell_selector);
      if (sqlen != std::nullopt)
        out.push_back(Long_edge_with_length{e, sqlen.value()});
    };

    bool collected = false;
    // C1: the fused edge pass has already applied `keep` to every finite edge,
    // in the same scan that answered resolution_reached() and collected the
    // collapse candidates. Three traversals of the edge set become one.
    if (m_precollected != nullptr)
    {
      long_edges_with_lengths.reserve(m_precollected->size());
      for (const auto& el : *m_precollected)
        long_edges_with_lengths.push_back(Long_edge_with_length{el.first, el.second});
      collected = true;
    }
#ifdef CGAL_LINKED_WITH_TBB
    if (!collected)
    if constexpr (std::is_convertible_v<typename Tr::Concurrency_tag, CGAL::Parallel_tag>)
    {
      long_edges_with_lengths
        = parallel_collect_finite_edges<Long_edge_with_length>(tr, keep);
      collected = true; // an empty result is a result, not a fallback
    }
#endif
    if (!collected)
      for (Edge e : tr.finite_edges())
        keep(e, long_edges_with_lengths);

    // longest first; stable to match the original bimap's ordering.
    //
    // `CGAL_TR_PARALLEL_SORT=1` (C6, replay of queue item R21) sorts in
    // parallel instead. tbb::parallel_sort is NOT stable, so ties have to be
    // broken explicitly or the comparator is not a strict weak ordering and
    // the sort is undefined. The tie-break is the owning cell's address: it is
    // a total order within a run, which is all a sort needs, and the candidate
    // ORDER ACROSS RUNS is already non-deterministic on the parallel path
    // because the collection is a parallel cell scan. Split's dependence is on
    // longest-FIRST, which both forms preserve exactly; R2 showed what
    // breaking that costs (-6.4%).
#ifdef CGAL_LINKED_WITH_TBB
    if (Parallel_tuning::get().parallel_candidate_sort)
    {
      tbb::parallel_sort(long_edges_with_lengths.begin(), long_edges_with_lengths.end(),
                         [](const Long_edge_with_length& a, const Long_edge_with_length& b) {
                           if (a.sqlength != b.sqlength)
                             return a.sqlength > b.sqlength;
                           return &*a.edge.first < &*b.edge.first;
                         });
    }
    else
#endif
    std::stable_sort(long_edges_with_lengths.begin(), long_edges_with_lengths.end(),
                     [](const Long_edge_with_length& a, const Long_edge_with_length& b) {
                       return a.sqlength > b.sqlength;
                     });

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    {
      std::ofstream ofs("long_edges.polylines.txt");
      for (const auto& le : long_edges_with_lengths)
        ofs << "2 " << point(le.edge.first->point())
            << " " << point(le.edge.second->point()) << std::endl;
    }
    m_can_be_split_ofs.open("can_be_split_edges.polylines.txt");
    m_split_failed_ofs.open("split_failed.polylines.txt");
    m_midpoints_ofs.open("midpoints.off");
    m_midpoints_ofs << "OFF" << std::endl;
    m_midpoints_ofs << long_edges_with_lengths.size() << " 0 0" << std::endl;
#endif

    Long_edges long_edges;
    long_edges.reserve(long_edges_with_lengths.size());
    for(const auto& ef : long_edges_with_lengths)
      long_edges.push_back(make_vertex_pair(ef.edge));
    return long_edges;
  }

  bool execute_operation(const Element_type& element, C3t3& c3t3) override
  {
    Tr& tr = c3t3.triangulation();
    const Edge_vv& e = element;

#ifdef CGAL_TR_MVLZ_PROBE
    // The window opens BEFORE is_edge(). It used to open after, which left
    // is_edge()'s walk over the whole star of the first endpoint outside every
    // measurement -- the one part of the operation whose safety had to be
    // argued instead of measured. The cost is that a stale candidate now
    // counts as an operation; it writes nothing, so it does not disturb the
    // write results.
    mvlz_reporter();
    Mvlz_probe<Tr> mvlz(tr);
    mvlz.begin("split", e.first, e.second);
    struct Mvlz_end {
      Mvlz_probe<Tr>& p; bool ok = false;
      ~Mvlz_end() { p.end(ok); }
    } mvlz_end{mvlz};
#endif

    Cell_handle cell;
    int i1, i2;
    if (!tr.tds().is_edge(e.first, e.second, cell, i1, i2))
      return false;

    Edge edge(cell, i1, i2);

    // check that splittability has not changed
    if (!can_be_split(edge, c3t3, m_protect_boundaries, m_cell_selector).can_be_split)
      return false;
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    m_can_be_split_ofs << "2 " << edge.first->vertex(edge.second)->point()
                        << " " << edge.first->vertex(edge.third)->point() << std::endl;
#endif

    m_visitor.before_split(tr, edge);
    Vertex_handle vh = split_edge(edge, m_cell_selector, c3t3);

    if (vh == Vertex_handle())
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      m_split_failed_ofs << "2 " << edge.first->vertex(edge.second)->point() << " "
                         << edge.first->vertex(edge.third)->point() << std::endl;
#endif
      return false;
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    m_midpoints_ofs << vh->point() << std::endl;
#endif
    m_visitor.after_split(tr, vh);
#ifdef CGAL_TR_MVLZ_PROBE
    mvlz_end.ok = true;
#endif
    return true;
  }

  /**
  * Splitting `element` rewrites the cells incident to its two vertices, so
  * both stars must be held. Only used by the parallel executor.
  */
  bool lock_zone(const Element_type& element, const C3t3& c3t3) const
  {
    const Tr& tr = c3t3.triangulation();

    const int mode = Parallel_tuning::get().mvlz_split_zone;
    if (mode == 2)
    {
      // SABOTAGE: the two endpoints only, which cannot possibly cover the
      // cells the split rewires. Exists so the lock-coverage check can be
      // shown able to report a violation on this binary.
      return tr.try_lock_vertex(element.first) && tr.try_lock_vertex(element.second);
    }
    if (mode == 1)
      return lock_zone_mvlz(element, tr);

    std::vector<Cell_handle> inc_cells_first, inc_cells_second;
    Star_census::hit(7, 0);
    return tr.try_lock_and_get_incident_cells(element.first, inc_cells_first)
        && tr.try_lock_and_get_incident_cells(element.second, inc_cells_second);
  }

  /**
  * The MEASURED zone (MVLZ_SPLIT.md): the vertices of the cells incident to
  * the edge, plus one apex per ring facet.
  *
  * Why this is the whole footprint. A Valgrind Lackey trace of 92.5M memory
  * accesses over 10 splits recorded every load and store the operation
  * performs. Every cell and vertex of the triangulation was in the trace's
  * identity manifest, so a miss would have been reported rather than lost:
  *  - cells incident to the edge      3103 loads,  136 stores
  *  - their facet-neighbours          1313 loads,  170 stores
  *  - EVERY OTHER CELL                   0 loads,    0 stores
  * and nothing at all beyond graph distance 1. An independent snapshot/diff
  * over 741 splits on two meshes agrees on the write half.
  *
  * Three things this relies on, none of them obvious:
  *
  * 1. `execute_operation()` calls `tds().is_edge()` BEFORE anything this
  *    measurement covered, and that walks the whole star of the first
  *    endpoint. It is safe without holding that star because we hold the
  *    endpoint itself, and a thread may only modify a cell while holding all
  *    four of its vertices -- so no cell incident to a vertex we hold can be
  *    modified under us. That is the protocol invariant, not a measurement,
  *    and it is the main thing the crash soak is testing.
  *
  * 2. The mirror cell across a ring facet shares three vertices with the ring
  *    cell, which are already held; only its apex is missing. Same argument as
  *    the shipped A1/A2 apex-only halo.
  *
  * 3. `is_edge()` failing means the candidate went stale, NOT that the zone is
  *    contended. It must return true here: the executor spins
  *    `while (!lock_zone(...))`, so returning false for a stale candidate is
  *    an infinite loop. execute_operation() re-runs is_edge() and drops it.
  */
  bool lock_zone_mvlz(const Element_type& element, const Tr& tr) const
  {
    if (!tr.try_lock_vertex(element.first) || !tr.try_lock_vertex(element.second))
      return false;

    // The split CREATES a vertex at the midpoint, and the lock grid is keyed
    // on position, so the midpoint's grid cell must be held too -- the same
    // reason collapse locks its destination (the 2026-08-31 root-cause crash).
    // The lock-coverage check shows split has never done this: the shipped
    // both-stars zone leaves the new vertex unheld in 180 of 5008 changed
    // cells, purely because the midpoint usually happens to fall in a grid
    // cell some star vertex already covers. That is luck, not protection, and
    // a smaller zone gets less of it.
    if (!tr.try_lock_point(CGAL::midpoint(point(element.first->point()),
                                          point(element.second->point()))))
      return false;

    // is_edge() is not a read. TDS::is_edge() marks tds_data() on every cell
    // of the FIRST endpoint's star and clears it on scope exit, so it writes
    // shared scratch state on the whole star. That write is invisible to a
    // snapshot/diff because it restores the old value, which is exactly why
    // the first version of this zone measured clean and was still wrong: the
    // Lackey trace caught 54 stores per 8 splits on star cells that are
    // neither ring nor mirror. Two threads resolving overlapping edges would
    // corrupt each other's marks. So the first endpoint's star must be held
    // in full; only the second endpoint's exclusive star is saved.
    std::vector<Cell_handle> star_first;
    Star_census::hit(7, 0);
    if (!tr.try_lock_and_get_incident_cells(element.first, star_first))
      return false;

    Cell_handle c;
    int i1, i2;
    if (!tr.tds().is_edge(element.first, element.second, c, i1, i2))
      return true;                       // stale, not contended -- see (3)

    typename Tr::Cell_circulator circ = tr.incident_cells(Edge(c, i1, i2));
    const typename Tr::Cell_circulator end = circ;
    do
    {
      const Cell_handle rc = circ;
      if (!tr.try_lock_cell(rc))         // the ring cell's four vertices
        return false;
      for (int k = 0; k < 4; ++k)
      {
        const Cell_handle n = rc->neighbor(k);
        if (!tr.try_lock_vertex(n->vertex(n->index(rc))))   // its apex only
          return false;
      }
    }
    while (++circ != end);

    return true;
  }

  // Where the element sits, for the spatial grouping the parallel executor
  // can do. One endpoint is enough: the two are one edge apart, far closer
  // than a grid cell.
  typename Tr::Geom_traits::Point_3 point_on_element(const Element_type& e) const
  {
    return typename Tr::Geom_traits().construct_point_3_object()(e.first->point());
  }

  // Opts out of lock elision: returning nothing means "always take the locks".
  void locked_vertices(const Element_type&,
                       boost::container::small_vector<Vertex_handle, 2>&) const {}

  void elision_vertices(const Element_type&,
                        boost::container::small_vector<Vertex_handle, 2>&) const {}

  static constexpr int zone_ring = 1;

  // longest edge first is the point of the ordering built in get_elements()
  static constexpr bool requires_ordered_processing = true;

  std::string operation_name() const override { return "Split long edges"; }

private:
  const std::vector<Edge_with_length>* m_precollected = nullptr;
};

} // internal
} // Tetrahedral_remeshing
} // CGAL

#endif // CGAL_INTERNAL_SPLIT_LONG_EDGES_H
