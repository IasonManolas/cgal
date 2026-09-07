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

#ifndef CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H
#define CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <chrono>

#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>
#include <utility>
#include <unordered_set>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/utility.h>
#include <CGAL/Tetrahedral_remeshing/internal/Elementary_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/MVLZ_probe.h>
#include <atomic>
#include <iostream>

#ifdef CGAL_LINKED_WITH_TBB
#include <boost/unordered/concurrent_flat_map.hpp>
#include <tbb/concurrent_priority_queue.h>
#include <tbb/concurrent_queue.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>
#include <thread>
#endif

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#include <CGAL/Real_timer.h>
#endif

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{


enum Edge_type     { FEATURE, BOUNDARY, INSIDE, MIXTE,
                     NO_COLLAPSE, INVALID, IMAGINARY, MIXTE_IMAGINARY, HULL_EDGE };
enum Collapse_type { TO_MIDPOINT, TO_V0, TO_V1, IMPOSSIBLE };
enum Result_type   { VALID,
                     V_PROBLEM, C_PROBLEM, E_PROBLEM,
                     ANGLE_PROBLEM,
                     TOPOLOGICAL_PROBLEM, ORIENTATION_PROBLEM, SHARED_NEIGHBOR_PROBLEM };

// A cell `c` incident to the edge being collapsed, with the two cells `n0`/`n1`
// that will take its place once it is removed, and the index of `c` in each of
// them. Collapsing the edge amounts to gluing `n0` and `n1` to each other
// across those indices.
template<typename Cell_handle>
struct Collapse_star_cell
{
  Cell_handle c;
  Cell_handle n0, n1;
  int c_in_n0, c_in_n1;

  // The two cells would end up pointing at the infinite vertex across from one
  // another, which is not a valid triangulation.
  template<typename Tr>
  bool has_infinite_adjacency(const Tr& tr) const
  {
    return tr.is_infinite(n0->vertex(c_in_n0))
        && tr.is_infinite(n1->vertex(c_in_n1));
  }
};

// `c` must be a cell incident to the edge (`v0`, `v1`) being collapsed, and
// this must be called before any neighbor around that edge is rewired:
// `index()` looks `c` up in the neighbor array of `n0`/`n1`, so it would no
// longer find it once `set_neighbor()` has run.
template<typename CellRef, typename Vertex_handle>
auto make_collapse_star_cell(CellRef c, Vertex_handle v0, Vertex_handle v1)
{
  using Cell_handle = std::decay_t<decltype(c->neighbor(0))>;

  const Cell_handle n0 = c->neighbor(c->index(v0));
  const Cell_handle n1 = c->neighbor(c->index(v1));

  return Collapse_star_cell<Cell_handle>{c, n0, n1, n0->index(c), n1->index(c)};
}

template<typename C3t3>
class CollapseTriangulation
{
  typedef typename C3t3::Triangulation                        Tr;
  typedef typename C3t3::Edge                                 Edge;
  typedef typename C3t3::Cell_handle                          Cell_handle;
  typedef typename C3t3::Vertex_handle                        Vertex_handle;
  typedef typename C3t3::Subdomain_index                      Subdomain_index;
  typedef typename C3t3::Surface_patch_index                  Surface_patch_index;
  typedef typename C3t3::Triangulation::Point                 Point_3;
  typedef typename C3t3::Triangulation::Geom_traits::Vector_3 Vector_3;

public:
  CollapseTriangulation(const Edge& e,
                        const std::unordered_set<Cell_handle>& cells_to_insert,
                        Collapse_type _collapse_type)
    : collapse_type(_collapse_type)
    , v0_init(e.first->vertex(e.second))
    , v1_init(e.first->vertex(e.third))
  {
    typedef std::array<int, 3> Facet;
    typedef std::array<int, 4> Tet;

    /*vertex of main tr - vertex of collapse tr*/
    // the star holds around twenty distinct vertices, so a linear scan over
    // handles compares cheaper than hashing them : hashing a compact-container
    // iterator goes through Time_stamper, and this runs on every attempt
    boost::container::small_vector<std::pair<Vertex_handle, int>, 32> v2i;
    const auto index_of = [&v2i](const Vertex_handle vh) -> int
    {
      for (const std::pair<Vertex_handle, int>& p : v2i)
        if (p.first == vh)
          return p.second;
      return -1;
    };

    std::vector<Point_3> points;
    std::vector<Tet> finite_cells;
    std::vector<int> subdomains;
    finite_cells.reserve(cells_to_insert.size());
    subdomains.reserve(cells_to_insert.size());

    for (Cell_handle ch : cells_to_insert)
    {
      Tet tet;
      for (int i = 0; i < 4; ++i)
      {
        const Vertex_handle vh = ch->vertex(i);
        int id = index_of(vh);
        if (id == -1)
        {
          id = static_cast<int>(points.size());
          v2i.emplace_back(vh, id);
          points.push_back(vh->point());
        }
        tet[i] = id;
      }
      finite_cells.push_back(tet);
      subdomains.push_back(ch->subdomain_index());
    }

    CGAL_expensive_assertion(index_of(v1_init) != -1);
    CGAL_expensive_assertion(index_of(v0_init) != -1);

    // finished
    std::vector<Vertex_handle> new_vertices;
    std::map<Facet, typename C3t3::Surface_patch_index> border_facets;
    if (CGAL::SMDS_3::build_triangulation_impl(
            triangulation, points, finite_cells, subdomains, border_facets,
            new_vertices, /*verbose*/ false,
            /*replace_domain_0*/ false,
            /*allow_non_manifold*/false))
    {
      CGAL_expensive_assertion(triangulation.tds().is_valid());
      CGAL_assertion(triangulation.infinite_vertex() == new_vertices[0]);

      // update()
      vh0 = new_vertices[index_of(v0_init) + 1];
      vh1 = new_vertices[index_of(v1_init) + 1];

      Cell_handle ch;
      int i0, i1;
      not_an_edge = true;
      CGAL_assertion(triangulation.tds().is_vertex(vh0));
      CGAL_assertion(triangulation.tds().is_vertex(vh1));
      if (triangulation.is_edge(vh0, vh1, ch, i0, i1))
      {
        edge = Edge(ch, i0, i1);
        not_an_edge = false;
      }
    }
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    else
      std::cout << "Warning : CollapseTriangulation is not valid!" << std::endl;
#endif
  }

  Result_type collapse()
  {
    if (not_an_edge)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "CollapseTriangulation::Not an edge..." << std::endl;
#endif
      return E_PROBLEM;
    }
    else
    {
      Vector_3 v0_new_pos = vec(vh0->point());

      if (collapse_type == TO_MIDPOINT){
        v0_new_pos = v0_new_pos + 0.5 * Vector_3(point(vh0->point()), point(vh1->point()));
      }
      else if (collapse_type == TO_V1){
        v0_new_pos = vec(point(vh1->point()));
      }

      std::unordered_set<Cell_handle> invalid_cells;

      typedef typename Tr::Cell_circulator Cell_circulator;
      Cell_circulator circ = triangulation.incident_cells(edge);
      Cell_circulator done = circ;

      std::vector<Cell_handle> cells_to_remove;

      //Update the vertex before removing it
      std::vector<Cell_handle> find_incident;
      triangulation.incident_cells(vh0, std::back_inserter(find_incident));

      std::vector<Cell_handle> cells_to_update;
      triangulation.incident_cells(vh1, std::back_inserter(cells_to_update));

      do
      {
        const auto sc = make_collapse_star_cell(circ, vh0, vh1);

        if (sc.n0->has_neighbor(sc.n1))
          return SHARED_NEIGHBOR_PROBLEM;

        //Update neighbors before removing cell
        sc.n0->set_neighbor(sc.c_in_n0, sc.n1);
        sc.n1->set_neighbor(sc.c_in_n1, sc.n0);

        Subdomain_index si_n0 = sc.n0->subdomain_index();
        Subdomain_index si_n1 = sc.n1->subdomain_index();
        Subdomain_index si = circ->subdomain_index();

        if (si_n0 != si && si_n1 != si)
          return TOPOLOGICAL_PROBLEM;

        if (sc.has_infinite_adjacency(triangulation))
          return TOPOLOGICAL_PROBLEM;

        if ( triangulation.is_infinite(sc.n0)
             && triangulation.is_infinite(sc.n1)
             && !triangulation.is_infinite(circ))
          return TOPOLOGICAL_PROBLEM;

        cells_to_remove.push_back(sc.c);

        invalid_cells.insert(sc.c);

      } while (++circ != done);

      // compute and keep worst angle
      Dihedral_angle_cosine curr_max_cos
        = (std::max)(max_cos_dihedral_angle_in_range(triangulation, cells_to_remove, false),
                     max_cos_dihedral_angle_in_range(triangulation, cells_to_update, false));


      vh0->set_point(Point_3(v0_new_pos.x(), v0_new_pos.y(), v0_new_pos.z()));
      vh1->set_point(Point_3(v0_new_pos.x(), v0_new_pos.y(), v0_new_pos.z()));

      Vertex_handle infinite_vertex = triangulation.infinite_vertex();

      bool v0_updated = false;
      for (const Cell_handle& ch : find_incident)
      {
        if (invalid_cells.find(ch) == invalid_cells.end()) //valid cell
        {
          if (triangulation.is_infinite(ch))
            infinite_vertex->set_cell(ch);
          else {
            vh0->set_cell(ch);
            v0_updated = true;
          }
        }
      }

      //Update the vertex before removing it
      for (Cell_handle ch : cells_to_update)
      {
        if (invalid_cells.find(ch) == invalid_cells.end()) //valid cell
        {
          ch->set_vertex(ch->index(vh1), vh0);

          if (triangulation.is_infinite(ch))
            infinite_vertex->set_cell(ch);
          else {
            if (!v0_updated) {
              vh0->set_cell(ch);
              v0_updated = true;
            }
          }
        }
      }

      if (!v0_updated){
        std::cout << "CollapseTriangulation::PB i cell not valid!!!" << std::endl;
        return V_PROBLEM;
      }
      triangulation.tds().delete_vertex(vh1);

      //Removing cells
      for (Cell_handle ch : cells_to_remove){
        triangulation.tds().delete_cell(ch);
      }

      // check validity of cells
      for (Cell_handle cit : triangulation.finite_cell_handles())
      {
        if (!is_well_oriented(triangulation, cit))
          return ORIENTATION_PROBLEM;
      }

      // check angles
      for (Cell_handle cit : triangulation.finite_cell_handles())
      {
        auto max_cos_after_collapse = max_cos_dihedral_angle(triangulation, cit, false);
        if (      curr_max_cos < max_cos_after_collapse  // angles decreased
         && acceptable_max_cos < max_cos_after_collapse) // && angles go below acceptable bound
          return ANGLE_PROBLEM;
      }

      //int si_nb_vh0 = nb_incident_subdomains(vh0, c3t3);
      //int si_nb_vh1 = nb_incident_subdomains(vh1, c3t3);
      //int vertices_subdomain_nb_vh0 = std::max(si_nb_vh0, si_nb_vh1);
      //bool is_on_hull_vh0 = is_on_convex_hull(vh0, c3t3) || is_on_convex_hull(vh1, c3t3);

      //if( is_valid_for_domains() )
      return VALID;

      // return TOPOLOGICAL_PROBLEM;
    }
  }

protected:
  Tr triangulation;

  const Collapse_type collapse_type;

  const Vertex_handle v0_init;
  const Vertex_handle v1_init;

  Vertex_handle vh0;
  Vertex_handle vh1;

  Edge edge;

  bool not_an_edge;

  const Dihedral_angle_cosine acceptable_max_cos{0.995}; // 0.995 cos <=> 5.7 degrees
};



template<typename C3t3, typename CellSelector>
Collapse_type get_collapse_type(const typename C3t3::Edge& edge,
                                const C3t3& c3t3,
                                CellSelector cell_selector)
{
  bool update_v0 = false;
  bool update_v1 = false;
  get_edge_info(edge, update_v0, update_v1, c3t3, cell_selector);

  if (update_v0 && update_v1) return TO_MIDPOINT;
  else if (update_v0)         return TO_V1;
  else if (update_v1)         return TO_V0;
  else                        return IMPOSSIBLE;
}

//template<typename C3t3>
//Edge_type get_edge_type(const typename C3t3::Edge& edge,
//                        const C3t3& c3t3)
//{
//  typedef typename C3t3::Vertex_handle Vertex_handle;
//  typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
//  typedef typename C3t3::Subdomain_index Subdomain_index;

//  const Vertex_handle & v0 = edge.first->vertex(edge.second);
//  const Vertex_handle & v1 = edge.first->vertex(edge.third);

//  const int dim0 = c3t3.in_dimension(v0);
//  const int dim1 = c3t3.in_dimension(v1);

//  const bool is_v0_on_hull = is_on_convex_hull(v0, c3t3);
//  const bool is_v1_on_hull = is_on_convex_hull(v1, c3t3);

//  if (c3t3.is_in_complex(edge))
//    return FEATURE;

//  else if (dim0 == 3 && dim1 == 3)
//    return INSIDE;

//  else if (dim0 == 2 && dim1 == 2)
//  {
//    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
//    Cell_circulator done = circ;

//    std::vector<Subdomain_index> indices;
//    do
//    {
//      Subdomain_index current_si = circ->subdomain_index();

//      if (std::find(indices.begin(), indices.end(), current_si) == indices.end()) {
//        indices.push_back(current_si);
//      }

//      Subdomain_index si_n0 = circ->neighbor(circ->index(v0))->subdomain_index();
//      Subdomain_index si_n1 = circ->neighbor(circ->index(v1))->subdomain_index();
//      if (si_n0 == si_n1 && si_n0 != current_si)
//        return NO_COLLAPSE;

//    } while (++circ != done);

//    const std::size_t nb_si_v0 = nb_incident_subdomains(v0, c3t3);
//    const std::size_t nb_si_v1 = nb_incident_subdomains(v1, c3t3);

//    if (indices.size() >= (std::min)(nb_si_v0, nb_si_v1)) {
//      return BOUNDARY;
//    }
//  }

//  //std::cerr << "ERROR : get_edge_type did not return anything valid!" << std::endl;
//  return NO_COLLAPSE;
//}

template<typename C3t3>
bool is_valid_collapse(const typename C3t3::Edge& edge,
                       const C3t3& c3t3)
{
  typedef typename C3t3::Vertex_handle Vertex_handle;
  typedef typename C3t3::Cell_handle   Cell_handle;
  typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

  Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
  Cell_circulator done = circ;
  do
  {
    int v0_id = circ->index(v0);
    int v1_id = circ->index(v1);

    Cell_handle n0_ch = circ->neighbor(v0_id);
    Cell_handle n1_ch = circ->neighbor(v1_id);

    if (n0_ch->has_vertex(v0)
        || n1_ch->has_vertex(v1)
        || n0_ch->has_neighbor(n1_ch))
    {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      if (c3t3.is_in_complex(edge))
        ++nb_invalid_collapse_short;
#endif
      return false;
    }
  }
  while (++circ != done);

  return true;
}

// The cells of `star` keep a positive orientation once `v_moved` is at
// `new_pos`. Cells that also have `v_other` disappear with the collapse.
template<typename C3t3, typename CellRange>
bool collapse_keeps_orientations(const CellRange& star,
                                 const typename C3t3::Vertex_handle v_moved,
                                 const typename C3t3::Vertex_handle v_other,
                                 const typename C3t3::Triangulation::Geom_traits::Point_3& new_pos)
{
  typedef typename C3t3::Triangulation::Geom_traits::Point_3 Point;

  for (const auto& ch : star)
  {
    if (ch->has_vertex(v_other))
      continue;

    std::array<Point, 4> pts = { point(ch->vertex(0)->point()),
                                 point(ch->vertex(1)->point()),
                                 point(ch->vertex(2)->point()),
                                 point(ch->vertex(3)->point()) };
    pts[ch->index(v_moved)] = new_pos;
    if (CGAL::orientation(pts[0], pts[1], pts[2], pts[3]) != CGAL::POSITIVE)
      return false;
  }
  return true;
}

template<typename C3t3>
bool is_valid_collapse(const typename C3t3::Edge& edge,
                       const Collapse_type& collapse_type,
                       const typename C3t3::Triangulation::Point& new_pos,
                       const C3t3& c3t3)
{
  typedef typename C3t3::Vertex_handle        Vertex_handle;
  typedef typename C3t3::Cell_handle          Cell_handle;
  typedef typename C3t3::Triangulation::Point Point;

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
  const bool in_cx = c3t3.is_in_complex(edge);
  if (in_cx)
  {
    if (collapse_type == TO_MIDPOINT)
      nb_test_midpoint++;
    else if (collapse_type == TO_V1)
      nb_test_v1++;
    else
      nb_test_v0++;
  }
#endif

  if (collapse_type == TO_V1 || collapse_type == TO_MIDPOINT)
  {
    std::vector<Cell_handle> cells_to_check;
    c3t3.triangulation().finite_incident_cells(v0,
        std::back_inserter(cells_to_check));

    for (const Cell_handle& ch : cells_to_check)
    {
      if (!ch->has_vertex(v1))
      {
        //check orientation
        std::array<Point, 4> pts = { ch->vertex(0)->point(),
                                       ch->vertex(1)->point(),
                                       ch->vertex(2)->point(),
                                       ch->vertex(3)->point()};
        pts[ch->index(v0)] = new_pos;
        if (CGAL::orientation(point(pts[0]), point(pts[1]), point(pts[2]), point(pts[3]))
            != CGAL::POSITIVE)
        {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
          if (in_cx)
          {
            if (collapse_type == TO_MIDPOINT)
              nb_orientation_midpoint++;
            else
              nb_orientation_v1++;
          }
#endif
          return false;
        }
      }
    }
  }
  if (collapse_type == TO_V0 || collapse_type == TO_MIDPOINT)
  {
    std::vector<Cell_handle> cells_to_check;
    c3t3.triangulation().finite_incident_cells(v1,
        std::back_inserter(cells_to_check));

    for (const Cell_handle& ch : cells_to_check)
    {
      if (!ch->has_vertex(v0))
      {
        //check orientation
        std::array<Point, 4> pts = { ch->vertex(0)->point(),
                                       ch->vertex(1)->point(),
                                       ch->vertex(2)->point(),
                                       ch->vertex(3)->point() };
        pts[ch->index(v1)] = new_pos;
        if (CGAL::orientation(point(pts[0]), point(pts[1]), point(pts[2]), point(pts[3]))
            != CGAL::POSITIVE)
        {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
          if (in_cx)
          {
            if (collapse_type == TO_MIDPOINT)
              nb_orientation_midpoint++;
            else
              nb_orientation_v0++;
          }
#endif
          return false;
        }
      }
    }
  }

  return is_valid_collapse(edge, c3t3);
}

template<typename Facet, typename Vh>
bool facet_has_edge(const Facet& f, const Vh v0, const Vh v1)
{
  std::array<std::array<int, 2>, 3> edges = {{ {{1,2}}, {{2,3}}, {{3,1}} }};

  for (int i = 0; i < 3; ++i)
  {
    const std::array<int, 2>& ei = edges[i];
    if ( f.first->vertex((f.second + ei[0]) % 4) == v0
      && f.first->vertex((f.second + ei[1]) % 4) == v1)
      return true;
    if ( f.first->vertex((f.second + ei[0]) % 4) == v1
      && f.first->vertex((f.second + ei[1]) % 4) == v0)
      return true;
  }
  return false;
}

template<typename C3t3, typename CellSelector>
bool collapse_preserves_surface_star(const typename C3t3::Edge& edge,
                                     const C3t3& c3t3,
                                     const typename C3t3::Triangulation::Point& new_pos,
                                     const CellSelector& cell_selector)
{
  typedef typename C3t3::Triangulation       Tr;
  typedef typename C3t3::Vertex_handle       Vertex_handle;
  typedef typename C3t3::Facet               Facet;
  typedef typename Tr::Geom_traits::Vector_3 Vector_3;
  typedef typename Tr::Geom_traits::Point_3  Point_3;

  const Tr& tr = c3t3.triangulation();

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);
  if (c3t3.in_dimension(v0) != 2 || c3t3.in_dimension(v1) != 2)
    return true;//other cases should not be treated here

  typename Tr::Geom_traits gt = c3t3.triangulation().geom_traits();
  typename Tr::Geom_traits::Construct_opposite_vector_3
    opp = gt.construct_opposite_vector_3_object();
  typename Tr::Geom_traits::Compute_scalar_product_3
    product = gt.compute_scalar_product_3_object();
  typename Tr::Geom_traits::Construct_normal_3
    normal = gt.construct_normal_3_object();

  std::unordered_set<Facet, boost::hash<Facet>> facets;
  tr.finite_incident_facets(v0, std::inserter(facets, facets.end()));
  tr.finite_incident_facets(v1, std::inserter(facets, facets.end()));

// note : checking a 2nd ring of facets does not change the result
//  std::unordered_set<Facet, boost::hash<Facet>> ring2;
//  for (const Facet& f : facets)
//  {
//    for (int i = 1; i < 4; ++i)
//    {
//      Vertex_handle vi = f.first->vertex((f.second + i) % 4);
//      tr.finite_incident_facets(vi, std::inserter(ring2, ring2.end()));
//    }
//  }
//  facets.insert(ring2.begin(), ring2.end());

  Vector_3 reference_normal = CGAL::NULL_VECTOR;
  //Point_3 reference_c;
  for (const Facet& f : facets)
  {
    if (!is_boundary(c3t3, f, cell_selector))
      continue;
    if (facet_has_edge(f, v0, v1))
      continue; //this facet will collapse if collapse happens

    std::array<Point_3, 3> pts = {{ point(f.first->vertex((f.second + 1) % 4)->point()),
                                    point(f.first->vertex((f.second + 2) % 4)->point()),
                                    point(f.first->vertex((f.second + 3) % 4)->point()) }};
    if(f.second % 2 == 0)
      std::swap(pts[0], pts[1]);

    Vector_3 n_before_collapse = normal(pts[0], pts[1], pts[2]);

    const Facet& mf = tr.mirror_facet(f);
    bool do_opp = false;
    if (  c3t3.triangulation().is_infinite(mf.first)
      ||  c3t3.subdomain_index(mf.first) < c3t3.subdomain_index(f.first))
    {
      n_before_collapse = opp(n_before_collapse);
      do_opp = true;
    }

    if (reference_normal == CGAL::NULL_VECTOR)
    {
      //reference_c = CGAL::centroid(pts[0], pts[1], pts[2]);
      reference_normal = n_before_collapse;
    }

    // check after move
    for (int i = 0; i < 3; ++i)
    {
      const Vertex_handle vi = f.first->vertex((f.second + i + 1) % 4);
      if (vi == v0 || vi == v1)
      {
        if (f.second % 2 == 0)
        {
          if(i == 0)      pts[1] = point(new_pos);
          else if(i == 1) pts[0] = point(new_pos);
          else            pts[2] = point(new_pos);
        }
        else
          pts[i] = point(new_pos);
        break;
      }
    }

    Vector_3 n_after_collapse = normal(pts[0], pts[1], pts[2]);
    if(do_opp)
      n_after_collapse = opp(n_after_collapse);

    const double dotref = product(reference_normal, n_after_collapse);
    if(dotref < 0)
      return false;
    const double dot = product(n_before_collapse, n_after_collapse);
    if(dot < 0)
      return false;

//    if (dot * dotref < 0)
//    {
//      std::cout << "\ncollapse edge : ";
//      std::cout << point(v0->point()) << " " << point(v1->point()) << std::endl;
//      const auto vs = c3t3.triangulation().vertices(f);
//      const std::array<Point_3, 3> ps = { {point(vs[0]->point()),
//                                           point(vs[1]->point()),
//                                           point(vs[2]->point())} };
//      std::cout << "facet : ";
//      std::cout << ps[0] << " " << ps[1] << " " << ps[2] << std::endl;
//      const Point_3 c = CGAL::centroid(ps[0], ps[1], ps[2]);
//      std::cout << "n_before_collapse ";
//      std::cout << c << " " << (c + n_before_collapse) << std::endl;
//      std::cout << "n_after_collapse  ";
//      std::cout << c << " " << (c + n_after_collapse) << std::endl;
//      std::cout << "reference_normal  ";
//      std::cout << reference_c << " " << (reference_c + reference_normal) << std::endl;
//      std::cout << std::endl;
//    }
//    if (dotref < 0 || dot < 0)
//      return false;
  }

  return true;
}

template<typename C3t3, typename Sizing, typename CellSelector>
bool are_edge_lengths_valid(const typename C3t3::Edge& edge,
                            const C3t3& c3t3,
                            const Collapse_type& collapse_type,
                            const typename C3t3::Triangulation::Point& new_pos,
                            const Sizing& sizing,
                            const CellSelector& cell_selector)
{
  //SqLengthMap::key_type is Vertex_handle
  //SqLengthMap::value_type is double
  typedef typename C3t3::Triangulation::Geom_traits::FT FT;
  typedef typename C3t3::Edge                           Edge;
  typedef typename C3t3::Vertex_handle                  Vertex_handle;

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

  std::vector<Edge> inc_edges;
  if (collapse_type == TO_V1 || collapse_type == TO_MIDPOINT)
    c3t3.triangulation().finite_incident_edges(v0,
      std::back_inserter(inc_edges));
  if (collapse_type == TO_V0 || collapse_type == TO_MIDPOINT)
    c3t3.triangulation().finite_incident_edges(v1,
      std::back_inserter(inc_edges));

  typename C3t3::Index new_index{};
  int new_dim = -1;
  FT sizing_at_new_pos = FT(0);
  if (collapse_type == TO_V0)
  {
    new_index = c3t3.index(v0);
    new_dim = c3t3.in_dimension(v0);
    sizing_at_new_pos = sizing_at_vertex(v0, sizing, c3t3, cell_selector);
  }
  else if (collapse_type == TO_V1)
  {
    new_index = c3t3.index(v1);
    new_dim = c3t3.in_dimension(v1);
    sizing_at_new_pos = sizing_at_vertex(v1, sizing, c3t3, cell_selector);
  }
  else if (collapse_type == TO_MIDPOINT)
  {
    new_index = max_dimension_index(v0, v1);
    new_dim = (std::max)(v0->in_dimension(), v1->in_dimension());
#ifdef CGAL_AVERAGE_SIZING_AFTER_COLLAPSE
    sizing_at_new_pos = sizing_at_midpoint(edge, new_dim, new_index, sizing, c3t3, cell_selector);
#else
    sizing_at_new_pos = sizing(point(new_pos), new_dim, new_index);
    if(new_dim < 3 && sizing_at_new_pos == FT(0))
      sizing_at_new_pos = max_sizing_in_incident_cells(edge, sizing, c3t3, cell_selector);
#endif
  }
  else
    CGAL_assertion(false);

  boost::container::flat_set<Vertex_handle,
    std::less<Vertex_handle>,
    boost::container::small_vector<Vertex_handle, 64> > examined;

  for (const Edge& ei : inc_edges)
  {
    if (is_outside(ei, c3t3, cell_selector))
      continue;

    Vertex_handle vh = ei.first->vertex(ei.second);
    if (vh == v0 || vh == v1)
      vh = ei.first->vertex(ei.third);
    if (vh == v0 || vh == v1)
      continue;

    // a vertex adjacent to both extremities is met once per star, and the test
    // below reads nothing but `vh` : running it again cannot change its outcome
    if (!examined.insert(vh).second)
      continue;

    const FT sqlen = CGAL::squared_distance(point(new_pos), point(vh->point()));
    const FT sizing_at_vh = sizing_at_vertex(vh, sizing, c3t3, cell_selector);
    const FT sqhigh
        = CGAL::square(FT(4) / FT(3)) * (std::max)(CGAL::square(sizing_at_vh),
                                                   CGAL::square(sizing_at_new_pos));
    if (sqlen > sqhigh)
      return false;

      //if (adaptive){
      //  if (is_boundary_edge(ei) || is_hull_edge(ei)){
      //    if (sqlen_i > split_length)
      //      return false;
      //  }
      //  else if (sqlen_i > 4.*getAimedLength(ei, aimed_length) / 3.){// && is_in_complex(ei)  ){
      //    return false;
      //  }
      //}
      //else {
      //
      //}
  }

  return true;
}

template<typename C3t3>
void merge_surface_patch_indices(const typename C3t3::Facet& f1,
                                 const typename C3t3::Facet& f2,
                                 C3t3& c3t3)
{
  const bool in_cx_f1 = c3t3.is_in_complex(f1);
  const bool in_cx_f2 = c3t3.is_in_complex(f2);

  if (in_cx_f1 && !in_cx_f2)
  {
    typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(f1);
    f2.first->set_surface_patch_index(f2.second, patch);
  }
  else if (in_cx_f2 && !in_cx_f1)
  {
    typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(f2);
    f1.first->set_surface_patch_index(f1.second, patch);
  }
  else if(in_cx_f1 && in_cx_f2)
  {
    CGAL_assertion(c3t3.surface_patch_index(f1) == c3t3.surface_patch_index(f2));

    typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(f2);
    c3t3.remove_from_complex(f2);
    f2.first->set_surface_patch_index(f2.second, patch);
  }
}

template<typename C3t3, typename CellSelector, typename ShortEdgesBimap>
typename C3t3::Vertex_handle
collapse(const typename C3t3::Cell_handle ch,
         const int to, const int from,
         CellSelector& cell_selector,
         C3t3& c3t3,
         ShortEdgesBimap& short_edges)
{
  typedef typename C3t3::Triangulation Tr;
  typedef typename C3t3::Vertex_handle Vertex_handle;
  typedef typename C3t3::Cell_handle   Cell_handle;
  typedef typename C3t3::Facet         Facet;
  typedef typename Tr::Cell_circulator Cell_circulator;

  Tr& tr = c3t3.triangulation();

  Vertex_handle vkept = ch->vertex(to);
  const Vertex_handle vdeleted = ch->vertex(from);

  //Update the vertex before removing it
  std::vector<Cell_handle> incident_to_vkept;
  tr.incident_cells(vkept, std::back_inserter(incident_to_vkept));

  std::vector<Cell_handle> incident_to_vdeleted;
  tr.incident_cells(vdeleted, std::back_inserter(incident_to_vdeleted));

  // Resolve the whole star first, without modifying anything: rejecting the
  // collapse once some neighbors have been rewired would leave the
  // triangulation half-collapsed.
  boost::container::small_vector<Collapse_star_cell<Cell_handle>, 30> incident_to_edge;
  Cell_circulator circ = tr.incident_cells(ch, to, from);
  Cell_circulator done = circ;
  do
  {
    const auto sc = make_collapse_star_cell(circ, vkept, vdeleted);

    if (sc.has_infinite_adjacency(tr))
      return Vertex_handle();

    incident_to_edge.push_back(sc);
  }
  while (++circ != done);

  for (const auto& sc : incident_to_edge)
  {
    for (int i = 0; i < 4; ++i)
    {
      const Vertex_handle vi = sc.c->vertex(i);
      if (vi != vkept && vi != vdeleted)
      {
        const Facet fi(sc.c, i);
        if (c3t3.is_in_complex(fi))
          c3t3.remove_from_complex(fi);
      }
    }
  }

  if(c3t3.is_in_complex(ch->vertex(from), ch->vertex(to)))
    c3t3.remove_from_complex(ch->vertex(from), ch->vertex(to));

  std::vector<Cell_handle> cells_to_remove;
  std::unordered_set<Cell_handle> invalid_cells;

  for(const auto& sc : incident_to_edge)
  {
    //Merge surface patch indices
    merge_surface_patch_indices(Facet(sc.n0, sc.c_in_n0),
                                Facet(sc.n1, sc.c_in_n1),
                                c3t3);

    //Update neighbors before removing cell
    CGAL_TR_PROBE_CELL_WRITE(tr, sc.n0, "set_neighbor n0");
    CGAL_TR_PROBE_CELL_WRITE(tr, sc.n1, "set_neighbor n1");
    sc.n0->set_neighbor(sc.c_in_n0, sc.n1);
    sc.n1->set_neighbor(sc.c_in_n1, sc.n0);

    //Update vertices cell pointer
    for (int i = 0; i < 3; i++)
    {
      int vid = Tr::vertex_triple_index(sc.c_in_n0, i);
      CGAL_TR_PROBE_VERTEX_WRITE(tr, sc.n0->vertex(vid), "set_cell n0.vertex");
      sc.n0->vertex(vid)->set_cell(sc.n0);
    }
    for (int i = 0; i < 3; i++)
    {
      int vid = Tr::vertex_triple_index(sc.c_in_n1, i);
      CGAL_TR_PROBE_VERTEX_WRITE(tr, sc.n1->vertex(vid), "set_cell n1.vertex");
      sc.n1->vertex(vid)->set_cell(sc.n1);
    }

    cells_to_remove.push_back(sc.c);
    invalid_cells.insert(sc.c);
  }

  const Vertex_handle infinite_vertex = tr.infinite_vertex();

  bool v0_updated = false;
  for (const Cell_handle& c : incident_to_vkept)
  {
    if (invalid_cells.find(c) == invalid_cells.end())//valid cell
    {
      if (tr.is_infinite(c))
      {
        CGAL_TR_PROBE_VERTEX_WRITE(tr, infinite_vertex, "set_cell infinite (vkept loop)");
        infinite_vertex->set_cell(c);
      }
      //else {
      CGAL_TR_PROBE_VERTEX_WRITE(tr, vkept, "set_cell vkept");
      vkept->set_cell(c);
      v0_updated = true;
      //}
    }
  }

  // update complex edges
  for (const Cell_handle& c : incident_to_vdeleted)
  {
    for (const auto& ei : cell_edges(c, tr))
    {
      remove_from_bimap(ei, short_edges);

      const Vertex_handle eiv0 = c->vertex(ei.second);
      const Vertex_handle eiv1 = c->vertex(ei.third);
      if (eiv1 == vdeleted && eiv0 != vkept) //replace eiv1 by vkept
      {
        if (c3t3.is_in_complex(eiv0, eiv1))
        {
          if (!c3t3.is_in_complex(eiv0, vkept))
            c3t3.add_to_complex(eiv0, vkept, c3t3.curve_index(eiv0, eiv1));
          c3t3.remove_from_complex(eiv0, eiv1);
        }
      }
      else if (eiv0 == vdeleted && eiv1 != vkept) //replace eiv0 by vkept
      {
        if (c3t3.is_in_complex(eiv0, eiv1))
        {
          if (!c3t3.is_in_complex(vkept, eiv1))
            c3t3.add_to_complex(vkept, eiv1, c3t3.curve_index(eiv0, eiv1));
          c3t3.remove_from_complex(eiv0, eiv1);
        }
      }
    }
  }

  //Update the vertex before removing it
  for (const Cell_handle& c : incident_to_vdeleted)
  {
    if (invalid_cells.find(c) == invalid_cells.end()) //valid cell
    {
      CGAL_TR_PROBE_CELL_WRITE(tr, c, "set_vertex vdeleted->vkept");
      c->set_vertex(c->index(vdeleted), vkept);

      if (tr.is_infinite(c))
      {
        CGAL_TR_PROBE_VERTEX_WRITE(tr, infinite_vertex, "set_cell infinite (vdeleted loop)");
        infinite_vertex->set_cell(c);
      }
      //else {
      if (!v0_updated) {
        CGAL_TR_PROBE_VERTEX_WRITE(tr, vkept, "set_cell vkept (vdeleted loop)");
        vkept->set_cell(c);
        v0_updated = true;
      }
      //}
    }
  }

  if (!v0_updated)
    std::cout << "PB i cell not valid!!!" << std::endl;

  // Delete vertex
  CGAL_TR_PROBE_VERTEX_WRITE(tr, vdeleted, "delete_vertex");
  c3t3.triangulation().tds().delete_vertex(vdeleted);

  // Delete cells
  for (Cell_handle cell_to_remove : cells_to_remove)
  {
    // remove cell
    CGAL_TR_PROBE_CELL_WRITE(tr, cell_to_remove, "delete_cell");
    treat_before_delete(cell_to_remove, cell_selector, c3t3);
    c3t3.triangulation().tds().delete_cell(cell_to_remove);
  }

  return vkept;
}


template<typename C3t3, typename CellSelector, typename ShortEdgesBimap>
typename C3t3::Vertex_handle collapse(const typename C3t3::Edge& edge,
                                      const Collapse_type& collapse_type,
                                      CellSelector& cell_selector,
                                      C3t3& c3t3,
                                      ShortEdgesBimap& short_edges)
{
  typedef typename C3t3::Vertex_handle Vertex_handle;
  typedef typename C3t3::Triangulation::Point Point_3;

  Vertex_handle vh0 = edge.first->vertex(edge.second);
  Vertex_handle vh1 = edge.first->vertex(edge.third);

  const int dim_vh0 = c3t3.in_dimension(vh0);
  const int dim_vh1 = c3t3.in_dimension(vh1);

  Vertex_handle vh = Vertex_handle();

  const Point_3 p0 = vh0->point();
  const Point_3 p1 = vh1->point();

  //Collapse at mid point
  if (collapse_type == TO_MIDPOINT)
  {
    Point_3 new_position(CGAL::midpoint(point(vh0->point()), point(vh1->point())));
    vh0->set_point(new_position);
    vh1->set_point(new_position);

    vh = collapse(edge.first, edge.second, edge.third, cell_selector, c3t3, short_edges);
  }
  else //Collapse at vertex
  {
    if (collapse_type == TO_V1)
    {
      vh0->set_point(p1);
      vh = collapse(edge.first, edge.third, edge.second, cell_selector, c3t3, short_edges);
    }
    else //Collapse at v0
    {
      if (collapse_type == TO_V0)
      {
        vh1->set_point(p0);
        vh = collapse(edge.first, edge.second, edge.third, cell_selector, c3t3, short_edges);
      }
      else
        CGAL_assertion(false);
    }
  }

  // collapse() rejects an infinite adjacency before it rewires anything, so
  // the star is still the one we found. The two points are not : they were
  // moved above, in the expectation of a collapse that did not happen.
  if (vh == Vertex_handle())
  {
    vh0->set_point(p0);
    vh1->set_point(p1);
    return vh;
  }

  c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
  return vh;
}

template<typename C3t3>
bool is_cells_set_manifold(const C3t3&,
    std::unordered_set<typename C3t3::Cell_handle>& cells)
{
  typedef typename C3t3::Cell_handle Cell_handle;
  typedef typename C3t3::Vertex_handle Vh;
  typedef std::array<Vh, 3> FV;
  typedef std::pair<Vh, Vh> EV;

  // A facet is shared by exactly two cells, so it bounds the set when its
  // neighbor is outside : the triangulation already answers that, and asking
  // it costs one lookup of a cell handle where counting the facets of the set
  // meant hashing a triple of vertex handles for every facet of every cell.
  std::unordered_map<EV, int, boost::hash<EV>> edges;
  edges.reserve(4 * cells.size());

  for (Cell_handle c : cells)
  {
    for (int i = 0; i < 4; ++i)
    {
      if (cells.find(c->neighbor(i)) != cells.end())
        continue; // shared with another cell of the set

      const FV fvi = make_vertex_array(c->vertex((i + 1) % 4),
        c->vertex((i + 2) % 4),
        c->vertex((i + 3) % 4));

      for (int k = 0; k < 3; ++k)
      {
        const EV evi = make_vertex_pair(fvi[k], fvi[(k + 1) % 3]);
        typename std::unordered_map<EV, int, boost::hash<EV>>::iterator eit = edges.find(evi);
        if (eit == edges.end())
          edges.insert(std::make_pair(evi, 1));
        else
          eit->second++;
      }
    }
  }

  for (const auto& evv : edges)
    if (evv.second != 2)
      return false;

  return true;
}

enum Angle_verdict { ANGLES_REJECTED, ANGLES_ACCEPTED, ANGLES_UNDECIDED };

template<typename C3t3, typename CellRange>
Angle_verdict collapse_keeps_angles_acceptable(const typename C3t3::Edge& edge,
                                      const C3t3& c3t3,
                                      const Collapse_type collapse_type,
                                      const CellRange& star)
{
  using Tr = typename C3t3::Triangulation;
  using Cell_handle = typename Tr::Cell_handle;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Point_3 = typename Tr::Point;
  using Vector_3 = typename Tr::Geom_traits::Vector_3;
  using Subdomain_index = typename C3t3::Subdomain_index;

  const Dihedral_angle_cosine acceptable_max_cos{0.995}; // 0.995 cos <=> 5.7 degrees

  const Tr& tr = c3t3.triangulation();
  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

  Vector_3 new_pos = vec(v0->point());
  if (collapse_type == TO_MIDPOINT)
    new_pos = new_pos + 0.5 * Vector_3(point(v0->point()), point(v1->point()));
  else if (collapse_type == TO_V1)
    new_pos = vec(point(v1->point()));
  const auto p_new = point(Point_3(new_pos.x(), new_pos.y(), new_pos.z()));

  boost::container::flat_set<Cell_handle,
    std::less<Cell_handle>,
    boost::container::small_vector<Cell_handle, 32> > cells_to_remove;

  typename Tr::Cell_circulator circ = tr.incident_cells(edge);
  const typename Tr::Cell_circulator done = circ;
  do { cells_to_remove.insert(circ); } while (++circ != done);

  boost::container::small_vector<Cell_handle, 64> cells_to_update;
  tr.incident_cells(v1, std::back_inserter(cells_to_update));

  const Dihedral_angle_cosine curr_max_cos
    = (std::max)(max_cos_dihedral_angle_in_range(tr, cells_to_remove, false),
                 max_cos_dihedral_angle_in_range(tr, cells_to_update, false));

  const double angle_margin = 1e-9;
  const double acceptable_sq = acceptable_max_cos.signed_square_value();
  const double curr_max_sq = curr_max_cos.signed_square_value();
  bool undecided = false;

  const auto& gt = tr.geom_traits();
  for (const Cell_handle c : star)
  {
    if (cells_to_remove.find(c) != cells_to_remove.end())
      continue;
    if (tr.is_infinite(c) || c->subdomain_index() == Subdomain_index())
      continue;

    auto p_at = [&](const int i)
    {
      const Vertex_handle v = c->vertex(i);
      return (v == v0 || v == v1) ? p_new : point(v->point());
    };
    const Dihedral_angle_cosine after
      = max_cos_dihedral_angle(p_at(0), p_at(1), p_at(2), p_at(3), gt);
    const double after_sq = after.signed_square_value();

    if (CGAL::abs(after_sq - curr_max_sq) < angle_margin
     || CGAL::abs(after_sq - acceptable_sq) < angle_margin)
    {
      undecided = true;
      continue;
    }
    if (curr_max_cos < after && acceptable_max_cos < after)
      return ANGLES_REJECTED;
  }
  return undecided ? ANGLES_UNDECIDED : ANGLES_ACCEPTED;
}

template<typename C3t3,
         typename Sizing,
         typename CellSelector,
         typename ShortEdgesBimap,
         typename Visitor>
typename C3t3::Vertex_handle collapse_edge(const typename C3t3::Edge& edge,
    C3t3& c3t3,
    const Sizing& sizing,
    const bool /* protect_boundaries */,
    CellSelector cell_selector,
    ShortEdgesBimap& short_edges,
    Visitor& )
{
  typedef typename C3t3::Triangulation   Tr;
  typedef typename Tr::Point             Point;
  typedef typename Tr::Vertex_handle     Vertex_handle;
  typedef typename Tr::Cell_handle       Cell_handle;

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

  Collapse_type collapse_type = get_collapse_type(edge, c3t3, cell_selector);

#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
  const bool in_cx = c3t3.is_in_complex(edge);
  if (in_cx && collapse_type == IMPOSSIBLE)
    nb_impossible++;
#endif

  if (collapse_type == IMPOSSIBLE)
    return Vertex_handle();

  Point new_pos;
  switch(collapse_type)
  {
  case TO_V0:
    new_pos = v0->point(); break;
  case TO_V1:
    new_pos = v1->point(); break;
  default:
    CGAL_assertion(collapse_type == TO_MIDPOINT);
    new_pos = Point(CGAL::midpoint(point(v0->point()), point(v1->point())));
  }

  // The ring test depends on neither the collapse type nor the new position,
  // so it settles all three attempts below at once - and it is the cheap one:
  // one cell circulation, against a star walk plus an orientation predicate
  // per cell of the star.
  if (!is_valid_collapse(edge, c3t3))
    return Vertex_handle();

  // Each attempt below walks the star of v0 and/or of v1, and so do the angle
  // and manifold tests further down. The mesh is not touched in between, so
  // each star is walked once here and handed to all of them.
  boost::container::small_vector<Cell_handle, 64> star_v0, star_v1;
  bool has_star_v0 = false;
  bool has_star_v1 = false;
  const auto star_of_v0 = [&]() -> const boost::container::small_vector<Cell_handle, 64>&
  {
    if (!has_star_v0)
    {
      c3t3.triangulation().finite_incident_cells(v0, std::back_inserter(star_v0));
      has_star_v0 = true;
    }
    return star_v0;
  };
  const auto star_of_v1 = [&]() -> const boost::container::small_vector<Cell_handle, 64>&
  {
    if (!has_star_v1)
    {
      c3t3.triangulation().finite_incident_cells(v1, std::back_inserter(star_v1));
      has_star_v1 = true;
    }
    return star_v1;
  };

  const auto orientations_ok = [&](const Collapse_type ct, const Point& pos)
  {
    if ((ct == TO_V1 || ct == TO_MIDPOINT)
        && !collapse_keeps_orientations<C3t3>(star_of_v0(), v0, v1, point(pos)))
      return false;
    if ((ct == TO_V0 || ct == TO_MIDPOINT)
        && !collapse_keeps_orientations<C3t3>(star_of_v1(), v1, v0, point(pos)))
      return false;
    return true;
  };

  if (!orientations_ok(collapse_type, new_pos))
  {
    if (collapse_type == TO_MIDPOINT)
    {
      // with TO_MIDPOINT, we are authorized to test TO_V0 and TO_V1
      if (orientations_ok(TO_V0, v0->point()))
      {
        collapse_type = TO_V0;
        new_pos = v0->point();
      }
      else if (orientations_ok(TO_V1, v1->point()))
      {
        collapse_type = TO_V1;
        new_pos = v1->point();
      }
      else
      {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
        if (in_cx)
          nb_invalid_collapse++;
#endif
        return Vertex_handle();
      }
    }
    else
    {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      if (in_cx)
        nb_invalid_collapse++;
#endif
      return Vertex_handle();
    }
  }

  if (are_edge_lengths_valid(edge, c3t3, collapse_type, new_pos, sizing, cell_selector)
    && collapse_preserves_surface_star(edge, c3t3, new_pos, cell_selector))
  {
    CGAL_expensive_assertion(c3t3.triangulation().tds().is_edge(
                       edge.first->vertex(edge.second),
                       edge.first->vertex(edge.third)));

    std::unordered_set<Cell_handle> cells_to_insert;
    for (const Cell_handle ch : star_of_v0())
      cells_to_insert.insert(ch);
    for (const Cell_handle ch : star_of_v1())
      cells_to_insert.insert(ch);

    // the angle test is the one that discards most candidates, and the cheaper
    // of the two : it walks the star once, where is_cells_set_manifold() walks
    // the star of each of its vertices
    const Angle_verdict angles
      = collapse_keeps_angles_acceptable(edge, c3t3, collapse_type, cells_to_insert);
    if(angles == ANGLES_REJECTED)
      return Vertex_handle();

    if(!is_cells_set_manifold(c3t3, cells_to_insert))
      return Vertex_handle();

    if(angles == ANGLES_UNDECIDED)
    {
      CollapseTriangulation<C3t3> local_tri(edge, cells_to_insert, collapse_type);
      if(local_tri.collapse() != VALID)
        return Vertex_handle();
    }

#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
    if (in_cx)
      nb_valid_collapse++;
#endif
    return collapse(edge, collapse_type, cell_selector, c3t3, short_edges);
  }
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
  else if (in_cx)
    nb_invalid_lengths++;
#endif
  return Vertex_handle();
}

template<typename C3T3>
using Vertex_patch_cache = std::unordered_map<
    typename C3T3::Vertex_handle,
    std::optional<typename C3T3::Surface_patch_index> >;

// surface_patch_index(v) walks v's whole incident-facet star. The initial
// scan over all finite edges in collapse_short_edges() reaches the same
// vertex once per incident edge, and the scan does not modify the mesh, so
// the first answer stays valid for the rest of it.
#ifdef CGAL_TR_DIMSTATS
// N9's own gate: "count re-queue can_be_collapsed calls and their cache hit
// rate BEFORE building". Skipped once already, which is why the cache measured
// as a no-op with no explanation.
struct Patch_call_stats
{
  std::atomic<std::size_t> edges{0}, guard_passed{0}, walks{0}, hits{0};
  // Could the stored vertex index replace the walk entirely? Only if it is
  // (a) of the right kind -- it is a bare index whose MEANING depends on
  // in_dimension() -- and (b) still correct. collapse_short_edges.h contains
  // no call to set_index() at all, so nothing updates it after a collapse.
  std::atomic<std::size_t> dim0{0}, dim1{0}, dim2{0}, dim_other{0};
  std::atomic<std::size_t> idx_match{0}, idx_mismatch{0}, idx_no_patch{0};
  ~Patch_call_stats()
  {
    const std::size_t e = edges.load(), g = guard_passed.load(),
                      w = walks.load(), h = hits.load();
    std::cerr << "[patchstats] requeue_edges=" << e
              << " guard_passed=" << g
              << " (" << (e ? 100.0 * double(g) / double(e) : 0.0) << "% of edges)"
              << "  star_walks=" << w << "  cache_hits=" << h
              << " (" << ((w + h) ? 100.0 * double(h) / double(w + h) : 0.0) << "% hit)"
              << std::endl;
    const std::size_t d0 = dim0.load(), d1 = dim1.load(), d2 = dim2.load(),
                      dx = dim_other.load();
    std::cerr << "[patchidx] guard-passing endpoints by in_dimension:"
              << " dim0=" << d0 << " dim1=" << d1 << " dim2=" << d2
              << " other=" << dx << std::endl;
    std::cerr << "[patchidx] dim-2 endpoints, stored index vs walked patch:"
              << " match=" << idx_match.load()
              << " MISMATCH=" << idx_mismatch.load()
              << " walk_returned_nothing=" << idx_no_patch.load() << std::endl;
  }
};
inline Patch_call_stats& patch_stats() { static Patch_call_stats s; return s; }
#  define CGAL_TR_PATCHSTAT(f) (++::CGAL::Tetrahedral_remeshing::internal::patch_stats().f)
#else
#  define CGAL_TR_PATCHSTAT(f) ((void)0)
#endif

template<typename C3T3>
const std::optional<typename C3T3::Surface_patch_index>&
cached_surface_patch_index(const typename C3T3::Vertex_handle v,
                           const C3T3& c3t3,
                           Vertex_patch_cache<C3T3>& cache)
{
  const auto it = cache.find(v);
  if (it != cache.end())
  { CGAL_TR_PATCHSTAT(hits); return it->second; }

  CGAL_TR_PATCHSTAT(walks);
  return cache.emplace(v, surface_patch_index(v, c3t3)).first->second;
}

template<typename C3T3, typename CellSelector>
auto can_be_collapsed(const typename C3T3::Edge& e,
                      const C3T3& c3t3,
                      const bool protect_boundaries,
                      CellSelector cell_selector,
                      Vertex_patch_cache<C3T3>* patch_cache = nullptr)
{
  struct Collapsible
  {
    bool can_be_collapsed;
    bool on_boundary;
  };

  const bool in_cx = c3t3.is_in_complex(e);
  if(in_cx && protect_boundaries)
    return Collapsible{false, true /*boundary*/};

  const bool boundary = is_boundary(c3t3, e, cell_selector);
  if(boundary && protect_boundaries)
    return Collapsible{false, boundary};

  if(!is_selected(e, c3t3.triangulation(), cell_selector))
    return Collapsible{false, boundary};

  if(!boundary && !in_cx)
  {
    const auto v0 = e.first->vertex(e.second);
    const auto v1 = e.first->vertex(e.third);

    if(v0->in_dimension() != 3 && v1->in_dimension() != 3)
    {
      CGAL_TR_PATCHSTAT(guard_passed);
      CGAL_TR_PATCHSTAT(guard_passed);   // one per endpoint walked below
#ifdef CGAL_TR_DIMSTATS
      for (const auto& vv : { v0, v1 })
      {
        switch (vv->in_dimension())
        {
          case 0: CGAL_TR_PATCHSTAT(dim0); break;
          case 1: CGAL_TR_PATCHSTAT(dim1); break;
          case 2: CGAL_TR_PATCHSTAT(dim2); break;
          default: CGAL_TR_PATCHSTAT(dim_other); break;
        }
        if (vv->in_dimension() == 2)
        {
          const auto walked = surface_patch_index(vv, c3t3);
          if (walked == std::nullopt) CGAL_TR_PATCHSTAT(idx_no_patch);
          else
          {
            const auto stored =
              Mesh_3::internal::get_index<typename C3T3::Surface_patch_index>(vv->index());
            if (stored == walked.value()) CGAL_TR_PATCHSTAT(idx_match);
            else                          CGAL_TR_PATCHSTAT(idx_mismatch);
          }
        }
      }
#endif
      if (!patch_cache) { CGAL_TR_PATCHSTAT(walks); CGAL_TR_PATCHSTAT(walks); }
      const auto patch_v0 = patch_cache
        ? cached_surface_patch_index(v0, c3t3, *patch_cache)
        : surface_patch_index(v0, c3t3);
      const auto patch_v1 = patch_cache
        ? cached_surface_patch_index(v1, c3t3, *patch_cache)
        : surface_patch_index(v1, c3t3);

      if(patch_v0 != std::nullopt && patch_v1 != std::nullopt && patch_v0 != patch_v1)
        return Collapsible{false, boundary};
    }
  }

//   if(!is_internal(e, c3t3, cell_selector))
  return Collapsible {true, boundary};
}

#ifdef CGAL_LINKED_WITH_TBB
/**
* The parallel collapse work list. `collapse_edge()` reports the edges it
* destroys through `remove_from_bimap()`; a collapse running in parallel
* records them here instead, so that a thread reaching one of those candidates
* later skips it rather than following a cell handle that has been recycled.
* Edges are keyed by their vertex pair, not by the cell handle they were
* found through.
*/
/**
* What the parallel collapse keeps instead of the sequential work list.
*
* A collapse merges two vertices, and the one that goes takes with it every
* candidate that named it. Recording the *vertices* rather than the edges is
* what makes the re-queue safe: an edge is destroyed and re-created as its
* neighbourhood changes, but a vertex that has been merged away never comes
* back, so this set only ever vetoes candidates that really are gone.
*/
template<typename VertexHandle>
class Deleted_vertices
{
  // Open addressed on purpose. A split-ordered linked list (which is what
  // tbb::concurrent_unordered_set is) has to walk N scattered nodes and free
  // them one at a time when it is cleared, on one thread, at the end of every
  // collapse phase. This frees one contiguous block instead. Boost 1.83 has
  // concurrent_flat_map but not concurrent_flat_set, hence the map to bool.
  boost::concurrent_flat_map<VertexHandle, bool, boost::hash<VertexHandle> > m_vertices;

public:
  void insert(VertexHandle v) { m_vertices.emplace(v, true); }
  bool contains(VertexHandle v) const { return m_vertices.contains(v); }
};

/**
* `collapse_edge()` reports the edges it destroys so that a work list can drop
* them. The parallel collapse has no work list to keep up to date -- it finds
* out what is gone from `Deleted_vertices` -- so the reports go nowhere.
*/
struct No_work_list {};

template<typename Edge>
void remove_from_bimap(const Edge&, No_work_list&) {}
#endif // CGAL_LINKED_WITH_TBB

// The short edges left to collapse, shortest first. Edges are compared by
// their vertex pair, but stored with their orientation : `collapse_edge()`
// reads it to decide which extremity survives, so an edge already in the map
// keeps the orientation it entered with, and only its length is updated.
// The element side is only ever searched, never walked in order, so it is
// hashed : keeping it sorted meant comparing pairs of vertex handles
// O(log n) times for every edge the last collapse touched. The priority
// side keeps its order, which is what decides what runs next.
template<typename C3t3>
using Short_edges_bimap = boost::bimap<
    boost::bimaps::unordered_set_of<typename C3t3::Triangulation::Edge,
                                    Hash_edges<typename C3t3::Triangulation::Edge>,
                                    Equal_edges<typename C3t3::Triangulation::Edge> >,
    boost::bimaps::multiset_of<typename C3t3::Triangulation::Geom_traits::FT,
                               std::less<typename C3t3::Triangulation::Geom_traits::FT> > >;

#ifdef CGAL_TR_DIMSTATS
inline void dimstats(bool interior)
{
  struct Counts
  {
    std::atomic<std::size_t> total{0}, interior{0};
    ~Counts()
    {
      const std::size_t t = total.load(), i = interior.load();
      std::cerr << "[dimstats] collapses=" << t << " both_endpoints_interior=" << i
                << " (" << (t ? 100.0 * double(i) / double(t) : 0.0) << "%)"
                << std::endl;
    }
  };
  static Counts c;
  ++c.total;
  if (interior) ++c.interior;
}
#endif

template<typename C3t3,
         typename SizingFunction,
         typename CellSelector,
         typename Visitor>
class Edge_collapse_operation
    : public Elementary_operation<C3t3,
                                 typename C3t3::Triangulation::Edge,
                                 Short_edges_bimap<C3t3> >
{
public:
  using Tr = typename C3t3::Triangulation;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Edge = typename Tr::Edge;
  using FT = typename Tr::Geom_traits::FT;

  using Cell_handle = typename Tr::Cell_handle;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using Short_edges = Short_edges_bimap<C3t3>;
  using Base_operation = Elementary_operation<C3t3, Edge, Short_edges>;
  using Element_type = typename Base_operation::Element_type;
  using Element_range = typename Base_operation::Element_range;

private:
  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  Visitor& m_visitor;

#ifdef CGAL_LINKED_WITH_TBB
  Deleted_vertices<Vertex_handle> m_deleted_vertices; // parallel path only
#endif

public:
  Edge_collapse_operation(const SizingFunction& sizing,
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

  Element_range get_elements(const C3t3& c3t3) const override
  {
    Short_edges short_edges;
    const Tr& tr = c3t3.triangulation();

    // C1: the fused pass already ran can_be_collapsed + is_too_short over
    // every finite edge, in the same traversal that served split and
    // resolution_reached(). The bimap is still filled serially -- it is the
    // work list and its order decides what runs next.
    if (m_precollected != nullptr)
    {
      for (const Edge_with_length& el : *m_precollected)
        short_edges.insert(typename Short_edges::value_type(el.first, el.second));
      return short_edges;
    }

#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (std::is_convertible_v<typename Tr::Concurrency_tag, CGAL::Parallel_tag>)
    {
      {
        // The patch cache is deliberately not shared here: it is a plain map,
        // and the threads sharing nothing is the point.
        using Edge_with_length = std::pair<Edge, FT>;
        const std::vector<Edge_with_length> found
          = parallel_collect_finite_edges<Edge_with_length>(
              tr,
              [&](const Edge& e, std::vector<Edge_with_length>& out)
              {
                auto [collapsible, boundary]
                  = can_be_collapsed(e, c3t3, m_protect_boundaries, m_cell_selector);
                if (!collapsible)
                  return;
                const auto sqlen = is_too_short(e, boundary, m_sizing, c3t3, m_cell_selector);
                if (sqlen != std::nullopt)
                  out.emplace_back(e, sqlen.value());
              });

        // The bimap is filled serially: it is the work list, and its order is
        // what decides what runs next.
        for (const Edge_with_length& el : found)
          short_edges.insert(typename Short_edges::value_type(el.first, el.second));
        return short_edges;
      }
    }
#endif

    Vertex_patch_cache<C3t3> patch_cache;
    for (const Edge& e : tr.finite_edges())
    {
      auto [collapsible, boundary]
        = can_be_collapsed(e, c3t3, m_protect_boundaries, m_cell_selector, &patch_cache);
      if (!collapsible)
        continue;

      const auto sqlen = is_too_short(e, boundary, m_sizing, c3t3, m_cell_selector);
      if (sqlen != std::nullopt)
        short_edges.insert(typename Short_edges::value_type(e, sqlen.value()));
    }
    return short_edges;
  }

  /**
  * Unused: the sequential executor calls the three-argument overload below,
  * and the parallel executor calls execute_operation_vv(). Kept because
  * `Elementary_operation` declares it.
  */
  bool execute_operation(const Element_type& edge, C3t3& c3t3) override
  {
    Short_edges no_short_edges; // no work list to keep up to date
    return execute_operation(edge, c3t3, no_short_edges);
  }

#ifdef CGAL_LINKED_WITH_TBB
  /**
  * The parallel path works from a snapshot of vertex pairs rather than of
  * `Tr::Edge`s. An `Edge` names its edge through a cell handle, and a collapse
  * running on another thread may already have destroyed that cell, so
  * resolving one would follow recycled storage. Vertices survive until a
  * collapse merges them away, and every merged-away vertex is recorded in
  * `m_deleted_vertices` -- so a pair neither of whose vertices is in that set
  * still names two live vertices.
  */
  bool lock_zone(const Edge_vv& e, const C3t3& c3t3) const
  {
    if (is_gone(e))
    {
      CGAL_TR_ZS(stage_gone);
      return true; // nothing to lock; execute_operation_vv() will skip it
    }

    const Tr& tr = c3t3.triangulation();

    // SABOTAGE arm. The two endpoints and the destination, and nothing else:
    // no star, no halo. It cannot possibly cover the cells a collapse rewires,
    // and it exists so that "0 violations" from the measured arm can be shown
    // to be a result rather than a check that never fires. The midpoint is
    // kept because without it the arm fails for a second, unrelated reason
    // (the destination lock) and the two causes become impossible to tell
    // apart in the report.
    if (Parallel_tuning::get().mvlz_collapse_zone == 2)
    {
      if (!tr.try_lock_vertex(e.first) || !tr.try_lock_vertex(e.second))
        return false;
      if (is_gone(e)) return true;
      const auto& vs = tr.tds().vertices();
      if (!vs.is_used(e.first) || !vs.is_used(e.second)) return true;
      return tr.try_lock_point(CGAL::midpoint(point(e.first->point()),
                                              point(e.second->point())));
    }

    if (!tr.try_lock_vertex(e.first) || !tr.try_lock_vertex(e.second))
    {
      CGAL_TR_ZS(stage_endpoint);
      return false;
    }

    // Re-checked now that both vertices are held: another thread may have
    // merged one of them away between the test above and the locks.
    if (is_gone(e))
      return true;

    // Both vertices are held now, so the only thing that can still have
    // changed is whether their slots are live: a collapse that merged one of
    // them away before we took the locks. is_used() reads the slot's tag in
    // O(1) and without a lock, where is_vertex() would take the container's
    // block-list mutex and scan every block. A freed slot means the vertex is
    // gone for good, so this is a skip, not a retry -- retrying would spin
    // forever on a vertex that is never coming back.
    const auto& vertices = tr.tds().vertices();
    if (!vertices.is_used(e.first) || !vertices.is_used(e.second))
      return true;

    // A collapse does not only rewrite the star: it first MOVES the two
    // vertices -- to their midpoint, or one onto the other -- and only then
    // collapses. The spatial lock is keyed on where a vertex is when it is
    // locked, so as soon as one moves into a different grid cell nobody holds
    // that cell, and every write to the vertex and to the cells naming it is
    // unsynchronized from there on. Locking the two endpoints is therefore not
    // enough; the position they may move to has to be held as well. Moving one
    // endpoint onto the other lands on a point already held, so the midpoint
    // is the only destination left to take.
    if (!tr.try_lock_point(CGAL::midpoint(point(e.first->point()),
                                          point(e.second->point()))))
    {
      CGAL_TR_ZS(stage_midpoint);
      return false;
    }

    std::vector<Cell_handle> inc_cells_0, inc_cells_1;
    bool* const tls = zone_tls(tr);
    if (!tr.try_lock_and_get_incident_cells(e.first, inc_cells_0, tls))
    {
      CGAL_TR_ZS(stage_star0);
      return false;
    }
    Star_census::hit(6, inc_cells_0.size());
    if (!tr.try_lock_and_get_incident_cells(e.second, inc_cells_1, tls))
    {
      CGAL_TR_ZS(stage_star1);
      return false;
    }

    // ---- the halo ------------------------------------------------------
    // The claim this was written on: collapsing re-stitches the region it
    // removes to the cells around it, and in doing so calls set_cell() on the
    // vertices of those neighbouring cells -- including the one vertex of each
    // that lies outside the stars, at graph distance 2.
    //
    // MEASURED (MVLZ_COLLAPSE.md), and the claim does not survive. A Valgrind
    // Lackey trace of every load and store of eight collapses, with every cell
    // and vertex of the triangulation in the identity manifest so a miss would
    // be reported rather than lost:
    //
    //     depth 2, cell   :  88 loads    0 stores   (one byte offset)
    //     depth 2, vertex :   0 loads    0 stores
    //
    // Not one access of any kind reaches a vertex at distance 2, so the halo
    // locks vertices the operation never touches. The depth-2 CELLS it does
    // read are safe without it: a reader of a cell needs only one of its four
    // vertices, since a writer would need all four, and each of those cells
    // shares three vertices with a star cell we already hold. The snapshot/
    // diff agrees on the write half over 1,200 collapses on three meshes --
    // write radius 1, and every changed vertex at depth 0 or 1.
    //
    // So mode 1 drops it. Mode 0 keeps it, and is the control.
    if (Parallel_tuning::get().mvlz_collapse_zone != 1)
    {
      if (!lock_zone_halo(tr, inc_cells_0, inc_cells_1,
                          Parallel_tuning::get().apex_only_collapse_halo))
      {
        CGAL_TR_ZS(stage_halo);
        return false;
      }
    }
#ifdef CGAL_TR_ZONE_STATS
    // how many distinct lock-grid cells this zone actually spans
    {
      boost::container::small_vector<int, 64> seen;
      for (const std::vector<Cell_handle>* cs : { &inc_cells_0, &inc_cells_1 })
        for (const Cell_handle c : *cs)
          for (int k = 0; k < 4; ++k)
          {
            const int gi = tr.lock_grid_index(c->vertex(k)->point());
            if (gi >= 0 && std::find(seen.begin(), seen.end(), gi) == seen.end())
              seen.push_back(gi);
          }
      zone_stat_extent(seen.size());
    }
#endif

    // Self-check: the probe's own premise. If these fire, the probe and the
    // lock are not talking about the same thing, and no other report from it
    // means anything.
    CGAL_TR_PROBE_VERTEX_WRITE(tr, e.first,  "SELFCHECK lock_zone e.first");
    CGAL_TR_PROBE_VERTEX_WRITE(tr, e.second, "SELFCHECK lock_zone e.second");

    return true;
  }

  bool is_gone(const Edge_vv& e) const
  {
    return m_deleted_vertices.contains(e.first)
        || m_deleted_vertices.contains(e.second);
  }

  /**
  * Collapses `e`, and reports through `new_candidates` the edges that have
  * become short because of it. This is what the sequential path does by
  * updating its bimap in place: without it the parallel collapse stops at the
  * edges that were short when the snapshot was taken and leaves the mesh
  * coarsened less far than the sequential one.
  *
  * The edges reported are incident to the vertex the collapse kept, so they
  * lie in the zone this thread holds while it reads them.
  */
  template<typename OutputIterator>
  bool execute_operation_vv(const Edge_vv& e, C3t3& c3t3, OutputIterator new_candidates)
  {
    if (is_gone(e))
      return false;

#ifdef CGAL_TR_DIMSTATS
    // How many collapses could take a SMALL zone under a two-zone design?
    //
    // The depth-2 tds_data marking (MVLZ_COLLAPSE.md §0) is suspected to come
    // from can_be_collapsed() -> surface_patch_index(u) in the re-queue tail,
    // which is called only when BOTH endpoints of a re-queued edge have
    // in_dimension() != 3. Every re-queued edge is (vkept, u), so if vkept is
    // interior the call cannot happen for any u.
    //
    // vkept is interior whenever both endpoints are: the two merge rules in
    // this file disagree (max at collapse_type==TO_MIDPOINT, min in
    // merge_vertices) but 3 and 3 give 3 under either. So this predicate is
    // decidable at LOCK time from the two vertices already held, in O(1) --
    // unlike the k-ring elision, whose classification cost as much as the work
    // it avoided.
    {
      const bool interior = (e.first->in_dimension() == 3)
                         && (e.second->in_dimension() == 3);
      dimstats(interior);
    }
#endif

#ifdef CGAL_TR_MVLZ_PROBE
    // The window opens as early as it CAN, which is here: after is_gone() and
    // before is_edge(). Trap 1 of MVLZ_METHOD says open before the first line,
    // because the split probe first opened after is_edge() and a footprint
    // that excludes part of the operation is not the operation's footprint --
    // is_edge() marks tds_data() over a whole star, which is what went
    // missing. is_edge() is therefore inside the window here.
    //
    // is_gone() cannot be: it is the guard that says these two vertex handles
    // still name live vertices. A collapse DESTROYS one of its endpoints, so
    // a stale candidate's handles dangle, and the probe's first act is to walk
    // incident_cells() from both of them. Opening before is_gone() segfaults
    // on the first stale candidate -- observed, not predicted. Nothing is lost
    // by opening after it: is_gone() reads m_deleted_vertices, a hash set on
    // the side, and touches no cell and no vertex. This is the one place where
    // "before the first line" and "over the whole mesh footprint" differ, and
    // the second is the property that matters.
    //
    // The window also covers the RE-QUEUE TAIL below (finite_incident_edges
    // + can_be_collapsed + is_too_short around the kept vertex), because the
    // executor holds this zone until execute_operation_vv() returns. That
    // tail is part of what the zone has to protect whether or not it is part
    // of "the collapse", and it reads a 1-ring of a vertex that has just
    // MOVED, so leaving it outside the window would measure a different
    // operation from the one that runs.
    mvlz_reporter();
    Mvlz_probe<Tr> mvlz(c3t3.triangulation());
    mvlz.zone_today_has_apex_halo();          // collapse's shipped zone: A1
    // The two-zone predicate, decided from the two endpoints alone, in O(1),
    // at a point where both are held. If both are interior then the kept
    // vertex is interior (max and min of 3,3 are both 3, and the two merge
    // rules in this file disagree only elsewhere), so can_be_collapsed() can
    // never call surface_patch_index() for any re-queued edge.
    mvlz.classify((e.first->in_dimension() == 3 && e.second->in_dimension() == 3) ? 1 : 0);
    mvlz.begin("collapse", e.first, e.second);
    struct Mvlz_end {
      Mvlz_probe<Tr>& p; bool ok = false;
      ~Mvlz_end() { p.end(ok); }
    } mvlz_end{mvlz};
#endif

    Cell_handle cell;
    int i0, i1;
    if (!c3t3.triangulation().tds().is_edge(e.first, e.second, cell, i0, i1))
      return false;

    No_work_list no_work_list;
    const Vertex_handle vkept = collapse_edge(Edge(cell, i0, i1), c3t3, m_sizing,
                                              m_protect_boundaries, m_cell_selector,
                                              no_work_list, m_visitor);
    if (vkept == Vertex_handle())
      return false;

    // the collapse merged the two endpoints; the one that is not kept is gone
    m_deleted_vertices.insert(vkept == e.first ? e.second : e.first);

    std::vector<Edge> incident;
    c3t3.triangulation().finite_incident_edges(vkept, std::back_inserter(incident));

    // N9 -- a patch cache for the re-queue. Every edge here is (vkept, u), so
    // can_be_collapsed() computes surface_patch_index(vkept) once per incident
    // edge, ~20 times, for the same answer. The cache lives for THIS operation
    // only: the mesh is not modified again before the loop ends, which is the
    // condition the cache's own comment states, and a cache that outlived the
    // operation would be answering with patches from a mesh that has since
    // been collapsed.
    Vertex_patch_cache<C3t3> requeue_cache;
    Vertex_patch_cache<C3t3>* const cache_ptr =
      Parallel_tuning::get().requeue_patch_cache ? &requeue_cache : nullptr;

    for (const Edge& ei : incident)
    {
      CGAL_TR_PATCHSTAT(edges);
      const auto [collapsible, boundary]
        = can_be_collapsed(ei, c3t3, m_protect_boundaries, m_cell_selector, cache_ptr);
      if (!collapsible)
        continue;

      const std::optional<FT> sqlen
        = is_too_short(ei, boundary, m_sizing, c3t3, m_cell_selector);
      if (sqlen != std::nullopt)
        *new_candidates++ = std::make_pair(sqlen.value(), make_vertex_pair(ei));
    }
#ifdef CGAL_TR_MVLZ_PROBE
    mvlz_end.ok = true;
#endif
    return true;
  }

  // shortest edge first is the point of the ordering the bimap keeps
  static constexpr bool requires_ordered_processing = true;

#endif // CGAL_LINKED_WITH_TBB

private:
  // set_precollected() and get_elements() use this unconditionally, so it
  // cannot live inside the TBB guard: without TBB the class did not compile
  // at all. Only the C1 fused-edge-pass path ever sets it.
  const std::vector<Edge_with_length>* m_precollected = nullptr;
public:

  /**
  * Collapses `edge`, and keeps `short_edges` up to date : `collapse_edge()`
  * removes from it the edges it destroys, and the edges incident to the
  * vertex it keeps are re-evaluated here, since their length has changed.
  */
  bool execute_operation(const Element_type& edge, C3t3& c3t3,
                         Short_edges& short_edges)
  {
    const Vertex_handle vh = collapse_edge(edge, c3t3, m_sizing, m_protect_boundaries,
                                           m_cell_selector, short_edges, m_visitor);
    if (vh == Vertex_handle())
      return false;

    std::vector<Edge> incident_short;
    c3t3.triangulation().finite_incident_edges(vh, std::back_inserter(incident_short));
    for (const Edge& eshort : incident_short)
    {
      const auto [collapsible, boundary]
        = can_be_collapsed(eshort, c3t3, m_protect_boundaries, m_cell_selector);

      // an edge that can no longer be collapsed leaves the work list, rather
      // than being taken out of it later and refused by collapse_edge()
      std::optional<FT> sqlen;
      if (collapsible)
        sqlen = is_too_short(eshort, boundary, m_sizing, c3t3, m_cell_selector);

      update_bimap(eshort, short_edges, sqlen);
    }
    return true;
  }

  std::string operation_name() const override { return "Collapse short edges"; }
};

/**
* Collapse is the only operation whose elements change as it runs : collapsing
* an edge shortens the ones around the vertex it keeps, and destroys others.
* Its elements are therefore taken from a work list that `execute_operation()`
* keeps up to date, shortest first, rather than from a list collected once.
*/
template<typename C3t3,
         typename SizingFunction,
         typename CellSelector,
         typename Visitor>
class Elementary_operation_execution_sequential<
        Edge_collapse_operation<C3t3, SizingFunction, CellSelector, Visitor> >
{
  using Operation = Edge_collapse_operation<C3t3, SizingFunction, CellSelector, Visitor>;
  using Short_edges = typename Operation::Short_edges;
  using Edge = typename Operation::Edge;

public:
  bool execute(Operation& op, C3t3& c3t3) const
  {
    Short_edges short_edges = op.get_elements(c3t3);
    if (short_edges.empty())
      return false;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::size_t nb_done = 0;
    CGAL::Real_timer timer;
    timer.start();
#endif
    while (!short_edges.empty())
    {
      // the edge with shortest length
      typename Short_edges::right_map::iterator eit = short_edges.right.begin();
      const Edge e = eit->second;
      short_edges.right.erase(eit);

      if (op.execute_operation(e, c3t3, short_edges))
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        ++nb_done;
#endif
      }
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    timer.stop();
    std::cout << op.operation_name() << ": " << nb_done << " done ("
              << timer.time() << " sec)." << std::endl;
#endif
    return true;
  }
};

#ifdef CGAL_LINKED_WITH_TBB
/**
* The parallel counterpart. The bimap work list cannot be shared between
* threads, so the candidates are snapshotted from it once, shortest first, and
* drained from a concurrent queue in that order. What the sequential path does
* by re-queueing, this one does by recording the destroyed edges in the
* operation and skipping them.
*/
template<typename C3t3,
         typename SizingFunction,
         typename CellSelector,
         typename Visitor>
class Elementary_operation_execution_parallel<
        Edge_collapse_operation<C3t3, SizingFunction, CellSelector, Visitor> >
{
  using Operation = Edge_collapse_operation<C3t3, SizingFunction, CellSelector, Visitor>;
  using Short_edges = typename Operation::Short_edges;
  using Edge_vv = typename Operation::Edge_vv;
  using FT = typename Operation::FT;

  // an edge and its squared length, so the queue can keep the shortest on top
  using Candidate = std::pair<FT, Edge_vv>;
  struct Longer_first
  {
    bool operator()(const Candidate& a, const Candidate& b) const
    { return a.first > b.first; }
  };

public:
  bool execute(Operation& op, C3t3& c3t3) const
  {
    const Short_edges short_edges = op.get_elements(c3t3);
    if (short_edges.empty())
      return false;

    // Shortest first, and it must stay so as the collapses themselves shorten
    // further edges, so this is a priority queue rather than the plain queue
    // the other ordered operation uses.
    tbb::concurrent_priority_queue<Candidate, Longer_first> queue;
    for (auto it = short_edges.left.begin(); it != short_edges.left.end(); ++it)
      queue.push(Candidate(it->second, make_vertex_pair(it->first)));

    const std::size_t nb_candidates = queue.size();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    CGAL::Real_timer timer;
    timer.start();
#endif

    tbb::parallel_for(0, tbb::this_task_arena::max_concurrency(),
                      [&](int)
                      {
                        std::vector<Candidate> requeued;
                        Candidate candidate;
                        while (queue.try_pop(candidate))
                        {
                          const Edge_vv& e = candidate.second;
#ifdef CGAL_TR_ZONE_STATS
                          CGAL_TR_ZS(zones_locked);
                          CGAL_TR_ZS(zone_attempts);
                          std::size_t zs_fail = 0;
                          const bool zs_time =
                            ((Zone_stats::get().zone_attempts.load() & 63u) == 0u);
                          const auto zs_t0 = std::chrono::steady_clock::now();
#endif
                          while (!op.lock_zone(e, c3t3))
                          {
#ifdef CGAL_TR_ZONE_STATS
                            CGAL_TR_ZS(zone_attempts);
                            CGAL_TR_ZS(yields);
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
                          requeued.clear();
                          if (!op.execute_operation_vv(e, c3t3, std::back_inserter(requeued)))
                            CGAL_TR_ZS(zones_wasted);
                          c3t3.triangulation().unlock_all_elements();

                          for (const Candidate& c : requeued)
                            queue.push(c);
                        }
                      });

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    timer.stop();
    std::cout << op.operation_name() << ": " << nb_candidates
              << " candidates (" << timer.time() << " sec, parallel)." << std::endl;
#endif
    return true;
  }
};

#endif // CGAL_LINKED_WITH_TBB

}
}
}

#endif // CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H
