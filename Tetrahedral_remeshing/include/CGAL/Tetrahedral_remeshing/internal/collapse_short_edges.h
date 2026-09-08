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

#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>

#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <unordered_set>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/utility.h>
#include <CGAL/Tetrahedral_remeshing/internal/Elementary_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

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

    std::unordered_set<Vertex_handle> vertices_to_insert;
    for (Cell_handle ch : cells_to_insert)
    {
      for(int i = 0; i < 4; ++i)
        vertices_to_insert.insert(ch->vertex(i));
    }

    CGAL_expensive_assertion(vertices_to_insert.end()
      != std::find(vertices_to_insert.begin(), vertices_to_insert.end(), v1_init));
    CGAL_expensive_assertion(vertices_to_insert.end()
      != std::find(vertices_to_insert.begin(), vertices_to_insert.end(), v0_init));

    std::unordered_map<Vertex_handle, int> v2i;/*vertex of main tr - vertex of collapse tr*/

    //To add the vertices only once
    std::vector<Point_3> points;
    int index = 0;
    for (Vertex_handle vh : vertices_to_insert)
    {
      if (v2i.find(vh) == v2i.end())
      {
        points.push_back(vh->point());
        v2i.insert(std::make_pair(vh, index++));
      }
    }

    std::vector<Tet> finite_cells;
    std::vector<int> subdomains;
    for (Cell_handle ch : cells_to_insert)
    {
      finite_cells.push_back( { v2i.at(ch->vertex(0)),
                                v2i.at(ch->vertex(1)),
                                v2i.at(ch->vertex(2)),
                                v2i.at(ch->vertex(3)) } );
      subdomains.push_back(ch->subdomain_index());
    }

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
      vh0 = new_vertices[v2i.at(v0_init) + 1];
      vh1 = new_vertices[v2i.at(v1_init) + 1];

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
        int v0_id = circ->index(vh0);
        int v1_id = circ->index(vh1);

        Cell_handle n0_ch = circ->neighbor(v0_id);
        Cell_handle n1_ch = circ->neighbor(v1_id);

        int ch_id_in_n0 = n0_ch->index(circ);
        int ch_id_in_n1 = n1_ch->index(circ);

        if (n0_ch->has_neighbor(n1_ch))
          return SHARED_NEIGHBOR_PROBLEM;

        //Update neighbors before removing cell
        n0_ch->set_neighbor(ch_id_in_n0, n1_ch);
        n1_ch->set_neighbor(ch_id_in_n1, n0_ch);

        Subdomain_index si_n0 = n0_ch->subdomain_index();
        Subdomain_index si_n1 = n1_ch->subdomain_index();
        Subdomain_index si = circ->subdomain_index();

        if (si_n0 != si && si_n1 != si)
          return TOPOLOGICAL_PROBLEM;

        if ( triangulation.is_infinite(n0_ch->vertex(ch_id_in_n0))
             && triangulation.is_infinite(n1_ch->vertex(ch_id_in_n1)))
          return TOPOLOGICAL_PROBLEM;

        if ( triangulation.is_infinite(n0_ch)
             && triangulation.is_infinite(n1_ch)
             && !triangulation.is_infinite(circ))
          return TOPOLOGICAL_PROBLEM;

        cells_to_remove.push_back(circ);

        invalid_cells.insert(circ);

      } while (++circ != done);

      // the dihedral angles were compared before this copy was built, by
      // collapse_keeps_angles_acceptable(), on the very same cells and points :
      // whatever is left here has already passed that test

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

  boost::container::small_vector<Cell_handle, 30> incident_to_edge;
  Cell_circulator circ = tr.incident_cells(ch, to, from);
  Cell_circulator done = circ;
  do
  {
    for (int i = 0; i < 4; ++i)
    {
      const Vertex_handle vi = circ->vertex(i);
      if (vi != vkept && vi != vdeleted)
      {
        const Facet fi(circ, i);
        if (c3t3.is_in_complex(fi))
          c3t3.remove_from_complex(fi);
      }
    }
    incident_to_edge.push_back(circ);
  }
  while (++circ != done);

  if(c3t3.is_in_complex(ch->vertex(from), ch->vertex(to)))
    c3t3.remove_from_complex(ch->vertex(from), ch->vertex(to));

  std::vector<Cell_handle> cells_to_remove;
  std::unordered_set<Cell_handle> invalid_cells;

  for(const Cell_handle& c : incident_to_edge)
  {
    const int v0_id = c->index(vkept);
    const int v1_id = c->index(vdeleted);

    Cell_handle n0_ch = c->neighbor(v0_id);
    Cell_handle n1_ch = c->neighbor(v1_id);

    const int ch_id_in_n0 = n0_ch->index(c);
    const int ch_id_in_n1 = n1_ch->index(c);

    //Merge surface patch indices
    merge_surface_patch_indices(Facet(n0_ch, ch_id_in_n0),
                                Facet(n1_ch, ch_id_in_n1),
                                c3t3);

    //Update neighbors before removing cell
    n0_ch->set_neighbor(ch_id_in_n0, n1_ch);
    n1_ch->set_neighbor(ch_id_in_n1, n0_ch);

    //Update vertices cell pointer
    for (int i = 0; i < 3; i++)
    {
      int vid = Tr::vertex_triple_index(ch_id_in_n0, i);
      n0_ch->vertex(vid)->set_cell(n0_ch);
    }
    for (int i = 0; i < 3; i++)
    {
      int vid = Tr::vertex_triple_index(ch_id_in_n1, i);
      n1_ch->vertex(vid)->set_cell(n1_ch);
    }

    if (tr.is_infinite(n0_ch->vertex(ch_id_in_n0))
      && tr.is_infinite(n1_ch->vertex(ch_id_in_n1)))
    {
      std::cout << "Collapse infinite issue!" << std::endl;
      return Vertex_handle();
    }
    cells_to_remove.push_back(c);
    invalid_cells.insert(c);
  }

  const Vertex_handle infinite_vertex = tr.infinite_vertex();

  bool v0_updated = false;
  for (const Cell_handle& c : incident_to_vkept)
  {
    if (invalid_cells.find(c) == invalid_cells.end())//valid cell
    {
      if (tr.is_infinite(c))
        infinite_vertex->set_cell(c);
      //else {
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
      remove_from_bimap(std::make_pair(ei.first->vertex(ei.second),
                                       ei.first->vertex(ei.third)),
                        short_edges);

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
      c->set_vertex(c->index(vdeleted), vkept);

      if (tr.is_infinite(c))
        infinite_vertex->set_cell(c);
      //else {
      if (!v0_updated) {
        vkept->set_cell(c);
        v0_updated = true;
      }
      //}
    }
  }

  if (!v0_updated)
    std::cout << "PB i cell not valid!!!" << std::endl;

  // Delete vertex
  c3t3.triangulation().tds().delete_vertex(vdeleted);

  // Delete cells
  for (Cell_handle cell_to_remove : cells_to_remove)
  {
    // remove cell
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
    c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
  }
  else //Collapse at vertex
  {
    if (collapse_type == TO_V1)
    {
      vh0->set_point(p1);
      vh = collapse(edge.first, edge.third, edge.second, cell_selector, c3t3, short_edges);
      c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
    }
    else //Collapse at v0
    {
      if (collapse_type == TO_V0)
      {
        vh1->set_point(p0);
        vh = collapse(edge.first, edge.second, edge.third, cell_selector, c3t3, short_edges);
        c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
      }
      else
        CGAL_assertion(false);
    }
  }
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
  // neighbour is outside : the triangulation already answers that, and asking
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

/**
* Does the collapse leave the dihedral angles of the star acceptable?
*
* `CollapseTriangulation` answers this only after building a local copy of the
* star and running the collapse on it, although the answer depends on the
* geometry alone : the cells that survive are the star minus the ring of the
* edge, with both extremities moved to the collapse point. Evaluating it here
* leaves that copy unbuilt whenever it would have been rejected - which is what
* happens to nearly half of the candidates that reach it.
*
* The comparison it performs is reproduced exactly, including the way the
* midpoint is computed, so that the two agree down to the last bit.
*/
template<typename C3t3, typename CellRange>
bool collapse_keeps_angles_acceptable(const typename C3t3::Edge& edge,
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

  // same expression as CollapseTriangulation::collapse()
  Vector_3 new_pos = vec(v0->point());
  if (collapse_type == TO_MIDPOINT)
    new_pos = new_pos + 0.5 * Vector_3(point(v0->point()), point(v1->point()));
  else if (collapse_type == TO_V1)
    new_pos = vec(point(v1->point()));
  const auto p_new = point(Point_3(new_pos.x(), new_pos.y(), new_pos.z()));

  boost::container::flat_set<Cell_handle,
    std::less<Cell_handle>,
    boost::container::small_vector<Cell_handle, 32> > ring;

  typename Tr::Cell_circulator circ = tr.incident_cells(edge);
  const typename Tr::Cell_circulator done = circ;
  do { ring.insert(circ); } while (++circ != done);

  // worst angle before : the ring, plus the star of the vertex that disappears
  Dihedral_angle_cosine curr_max_cos = max_cos_dihedral_angle_in_range(tr, ring, false);

  boost::container::small_vector<Cell_handle, 64> star_v1;
  tr.finite_incident_cells(v1, std::back_inserter(star_v1));
  const Dihedral_angle_cosine cos_v1
    = max_cos_dihedral_angle_in_range(tr, star_v1, false);
  if (curr_max_cos < cos_v1)
    curr_max_cos = cos_v1;

  // worst angle after : the cells of the star that the collapse keeps
  const auto& gt = tr.geom_traits();
  for (const Cell_handle c : star)
  {
    if (ring.find(c) != ring.end())
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

    if (curr_max_cos < after && acceptable_max_cos < after)
      return false;
  }
  return true;
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

  if (!is_valid_collapse(edge, collapse_type, new_pos, c3t3))
  {
    if (collapse_type == TO_MIDPOINT)
    {
      // with TO_MIDPOINT, we are authorized to test TO_V0 and TO_V1
      if (is_valid_collapse(edge, TO_V0, v0->point(), c3t3))
      {
        collapse_type = TO_V0;
        new_pos = v0->point();
      }
      else if (is_valid_collapse(edge, TO_V1, v1->point(), c3t3))
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

    Vertex_handle v0_init = edge.first->vertex(edge.second);
    Vertex_handle v1_init = edge.first->vertex(edge.third);

    std::unordered_set<Cell_handle> cells_to_insert;
    c3t3.triangulation().finite_incident_cells(v0_init,
      std::inserter(cells_to_insert, cells_to_insert.end()));
    c3t3.triangulation().finite_incident_cells(v1_init,
      std::inserter(cells_to_insert, cells_to_insert.end()));

    // the angle test is the one that discards most candidates, and the cheaper
    // of the two : it walks the star once, where is_cells_set_manifold() walks
    // the star of each of its vertices
    if(!collapse_keeps_angles_acceptable(edge, c3t3, collapse_type, cells_to_insert))
      return Vertex_handle();

    if(!is_cells_set_manifold(c3t3, cells_to_insert))
      return Vertex_handle();

    CollapseTriangulation<C3t3> local_tri(edge, cells_to_insert, collapse_type);

    Result_type res = local_tri.collapse();
    if (res == VALID)
    {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      if (in_cx)
        nb_valid_collapse++;
#endif
      return collapse(edge, collapse_type, cell_selector, c3t3, short_edges);
    }
  }
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
  else if (in_cx)
    nb_invalid_lengths++;
#endif
  return Vertex_handle();
}

template<typename C3T3, typename CellSelector>
auto can_be_collapsed(const typename C3T3::Edge& e,
                      const C3T3& c3t3,
                      const bool protect_boundaries,
                      CellSelector cell_selector)
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
      const auto patch_v0 = surface_patch_index(v0, c3t3);
      const auto patch_v1 = surface_patch_index(v1, c3t3);

      if(patch_v0 != std::nullopt && patch_v1 != std::nullopt && patch_v0 != patch_v1)
        return Collapsible{false, boundary};
    }
  }

//   if(!is_internal(e, c3t3, cell_selector))
  return Collapsible {true, boundary};
}

// The short edges left to collapse, shortest first. Edges are keyed by their
// vertex pair and compared regardless of orientation, but stored with their
// orientation : `collapse_edge()` reads it to decide which extremity
// survives, so an edge already in the map keeps the orientation it entered
// with, and only its length is updated.
//
// The key is a vertex pair rather than an Edge (Cell_handle, i, j) because
// collapsing destroys and recycles cells: an Edge key would have to be
// compared -- and so dereferenced -- long after the cell it names is gone.
// Vertex handles stay meaningful, and are re-resolved to a current Edge with
// tds().is_edge() at the point of use.
template<typename C3t3>
using Short_edges_bimap = boost::bimap<
    boost::bimaps::set_of<std::pair<typename C3t3::Triangulation::Vertex_handle,
                                    typename C3t3::Triangulation::Vertex_handle>,
                          Compare_vertex_pairs<typename C3t3::Triangulation::Vertex_handle> >,
    boost::bimaps::multiset_of<typename C3t3::Triangulation::Geom_traits::FT,
                               std::less<typename C3t3::Triangulation::Geom_traits::FT> > >;

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
  using Cell_handle = typename Tr::Cell_handle;
  using Edge = typename Tr::Edge;
  using FT = typename Tr::Geom_traits::FT;

  using Short_edges = Short_edges_bimap<C3t3>;
  using Base_operation = Elementary_operation<C3t3, Edge, Short_edges>;
  using Element_type = typename Base_operation::Element_type;
  using Element_range = typename Base_operation::Element_range;

private:
  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  Visitor& m_visitor;

public:
  Edge_collapse_operation(const SizingFunction& sizing,
                        const CellSelector& cell_selector,
                        const bool protect_boundaries,
                        Visitor& visitor)
      : m_sizing(sizing)
      , m_cell_selector(cell_selector)
      , m_protect_boundaries(protect_boundaries)
      , m_visitor(visitor) {}

  Element_range get_elements(const C3t3& c3t3) const override
  {
    Short_edges short_edges;
    const Tr& tr = c3t3.triangulation();

    auto eval = [&](const Edge& e, std::vector<std::pair<Edge, FT>>& out)
    {
      auto [collapsible, boundary]
        = can_be_collapsed(e, c3t3, m_protect_boundaries, m_cell_selector);
      if (!collapsible)
        return;

      const auto sqlen = is_too_short(e, boundary, m_sizing, c3t3, m_cell_selector);
      if (sqlen != std::nullopt)
        out.emplace_back(e, sqlen.value());
    };

#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
    // Parallel cell-scan candidate collection (see parallel_collect_finite_edges).
    // The bimap itself is still built serially below: boost::bimap is not
    // concurrency-safe for insertion, and its internal ordering doesn't depend
    // on collection order anyway (it sorts by length on the right map).
    std::vector<std::pair<Edge, FT>> short_edges_vec =
      parallel_collect_finite_edges<std::pair<Edge, FT>>(tr, eval);
#else
    std::vector<std::pair<Edge, FT>> short_edges_vec;
    for (const Edge& e : tr.finite_edges())
      eval(e, short_edges_vec);
#endif

    for (const auto& es : short_edges_vec)
      short_edges.insert(typename Short_edges::value_type(
        std::make_pair(es.first.first->vertex(es.first.second),
                       es.first.first->vertex(es.first.third)),
        es.second));

    return short_edges;
  }

  /**
  * Link vertices `vi` (`vi != v0 && vi != v1`) of every boundary-incident
  * facet of `edge`. This is exactly the vertex set `get_collapse_type()` ->
  * `topology_test()` reads via `nb_incident_subdomains(vi, ...)`
  * (`tetrahedral_remeshing_helpers.h`), which walks `vi`'s full
  * `incident_cells` star -- outside the endpoint-only claim the parallel wave
  * selector uses. Needed so the wave selector (in
  * `Elementary_operation_execution_parallel<Edge_collapse_operation>`) can
  * widen its claim past `v0`/`v1` to cover every star this edge's own
  * `collapse_edge()` call will actually read during the parallel phase.
  * Usually empty: interior edges (no incident boundary facet) never reach
  * `has_several_subdomains()` inside `topology_test()`.
  */
  boost::container::small_vector<Vertex_handle, 8>
  topology_test_link_vertices(const Element_type& edge, const C3t3& c3t3) const
  {
    boost::container::small_vector<Vertex_handle, 8> out;
    const Vertex_handle v0 = edge.first->vertex(edge.second);
    const Vertex_handle v1 = edge.first->vertex(edge.third);
    const auto& tr = c3t3.triangulation();

    auto fcirc = tr.incident_facets(edge);
    const auto fdone = fcirc;
    do
    {
      if (tr.is_infinite(fcirc->first))
        continue;
      if (is_boundary(c3t3, *fcirc, m_cell_selector))
      {
        for (int i = 1; i < 4; ++i)
        {
          const Vertex_handle vi = fcirc->first->vertex((fcirc->second + i) % 4);
          if (vi != v0 && vi != v1)
            out.push_back(vi);
        }
      }
    } while (++fcirc != fdone);

    return out;
  }

  bool lock_zone(const Element_type& edge, const C3t3& c3t3) const override
  {
    auto& tr = c3t3.triangulation();
    const Vertex_handle v0 = edge.first->vertex(edge.second);
    const Vertex_handle v1 = edge.first->vertex(edge.third);
    if (!(tr.try_lock_vertex(v0) && tr.try_lock_vertex(v1)))
      return false;
    const auto& tds = tr.tds();
    if (!tds.vertices().is_used(v0) || !tds.vertices().is_used(v1))
      return false;
    boost::container::small_vector<Cell_handle, 64> inc_cells_0, inc_cells_1;
    return tr.try_lock_and_get_incident_cells(v0, inc_cells_0)
        && tr.try_lock_and_get_incident_cells(v1, inc_cells_1);
  }

  bool requires_ordered_processing() const override { return true; }

  typename Tr::Geom_traits::Point_3 point_on_element(const Element_type& e) const
  {
    auto cp = typename Tr::Geom_traits().construct_point_3_object();
    return cp(e.first->vertex(e.second)->point());
  }

  bool execute_operation(const Element_type& edge, C3t3& c3t3) override
  {
    Short_edges no_short_edges; // no work list to keep up to date
    return execute_operation(edge, c3t3, no_short_edges);
  }

  /**
  * Collapses `edge`, and keeps `short_edges` up to date : `collapse_edge()`
  * removes from it the edges it destroys, and the edges incident to the
  * vertex it keeps are re-evaluated here, since their length has changed.
  *
  * `short_edges` is templated so the parallel executor can pass a
  * Short_edges_delta collector instead of the bimap itself, and replay the
  * recorded updates once its wave has joined (see
  * tetrahedral_remeshing_helpers.h).
  */
  /**
  * Re-evaluates every edge incident to `vh` and writes the result to the work
  * list. Runs `can_be_collapsed()` / `is_too_short()` on each edge (vh, x), and
  * BOTH of those walk the incident cells of BOTH endpoints of the edge they are
  * given -- so this traverses the star of every neighbour x of vh, well beyond
  * vh's own star. That is why it cannot run inside a parallel wave under the
  * cheap endpoint claim, and why the executor defers it (see
  * Short_edges_delta::refresh_vertices). Safe to call once the wave has joined,
  * and safe sequentially, where it is called inline from execute_operation().
  */
  template<typename ShortEdges>
  void refresh_incident_edges(const Vertex_handle vh, C3t3& c3t3,
                              ShortEdges& short_edges)
  {
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

      auto key = std::make_pair(eshort.first->vertex(eshort.second),
                                eshort.first->vertex(eshort.third));
      update_bimap(key, short_edges, sqlen);
    }
  }

  template<typename ShortEdgesOrDelta>
  bool execute_operation(const Element_type& edge, C3t3& c3t3,
                         ShortEdgesOrDelta& short_edges)
  {
    const Vertex_handle vh = collapse_edge(edge, c3t3, m_sizing, m_protect_boundaries,
                                           m_cell_selector, short_edges, m_visitor);
    if (vh == Vertex_handle())
      return false;

    // Under the parallel executor this must NOT run here : it walks the stars of
    // vh's neighbours (see refresh_incident_edges), which lie outside the claim
    // this wave member holds. Record the vertex and let the executor re-evaluate
    // once the wave has joined and the mesh is stable again.
    if constexpr (Defers_incident_refresh<ShortEdgesOrDelta>::value)
      short_edges.refresh_vertices.push_back(vh);
    else
      refresh_incident_edges(vh, c3t3, short_edges);

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
  using Cell_handle = typename C3t3::Cell_handle;

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
      const auto vpair = eit->second;
      short_edges.right.erase(eit);

      // the work list is keyed by vertex pair; resolve it to a current Edge,
      // keeping the orientation the pair was stored with
      Cell_handle cell;
      int i0, i1;
      if (!c3t3.triangulation().tds().is_edge(vpair.first, vpair.second, cell, i0, i1))
        continue;
      const Edge e(cell, i0, i1);

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

#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB

/**
* Collapse's work list is mutated as edges are collapsed: collapsing an edge
* destroys some of its neighbours and shortens others, which must be
* reflected in the shared shortest-first bimap before the next pop -- that
* pop-and-update step is inherently sequential. But the per-edge work itself
* (geometry validation, and the local triangulation edit) does not depend on
* other, non-overlapping edges, so it can run concurrently once a batch of
* "next shortest" candidates has been popped.
*
* Elements are processed in waves of up to COLLAPSE_WAVE_FACTOR * nthreads,
* selected shortest-first from the work list. Unlike split, the wave is not
* protected by lock_zone/retry: short edges cluster spatially, so a wave of
* the globally shortest ones puts every thread on the same small region,
* almost every lock attempt fails, and the retry loop degenerates into a
* spin storm -- orders of magnitude slower than sequential collapse.
*
* Instead the wave is chosen to be conflict-free: an edge joins it only if
* neither endpoint has been claimed by a wave member, where selecting an
* edge claims every vertex of its incident cells. Two wave members then
* share no cell, so their collapses read and write disjoint regions and need
* no locking whatsoever. Edges rejected for conflict stay in the work list
* and are reconsidered by the next wave, so shortest-first priority
* survives.
*
* Two things make this pay, and both are easy to get wrong:
* - Only the endpoints are tested, never the whole star. Testing the star of
*   every candidate costs more than the geometry it protects (~4.8 s to
*   expose ~1.4 s of work, on bear.mesh at 8 threads); testing endpoints is
*   equally safe, by the argument at the test itself, and confines the star
*   walk to accepted edges.
* - Waves must be large. Only a fraction of the candidates scanned survive
*   the conflict test, so a wave sized like split's leaves most threads idle
*   -- hence COLLAPSE_WAVE_FACTOR, well above ORDERED_WAVE_FACTOR.
*
* On bear.mesh (3 iterations, 8 threads) collapse drops from 36.5 s to
* 28.4 s against the sequential executor, and end-to-end remeshing from
* 174 s to 143 s (medians of repeated runs; single runs on this workload
* scatter by several seconds and one outlier is enough to invert the
* comparison, so measure it repeatedly). CGAL_TET_REMESHING_COLLAPSE_SEQUENTIAL
* opts back out.
*
* What is left is the serial replay below, still ~8.8 s of the last
* iteration -- the next thing worth attacking.
*
* Work-list updates are *not* applied to the bimap from the workers. A
* single collapse triggers on the order of a hundred of them (all six edges
* of every cell incident to the deleted vertex, plus every edge incident to
* the kept one), so pushing them through a shared lock costs more than the
* parallelism they sit inside -- measurably slower than running collapse
* sequentially. Each worker records its updates into a thread-local
* Short_edges_delta instead, and the executor replays them into the bimap
* after the wave joins, in per-thread order (see
* tetrahedral_remeshing_helpers.h).
*
* Wave members keep their Edge (Cell_handle, i, j) as stored in the work
* list -- including its orientation, which collapse_edge() reads to decide
* which extremity survives. That is safe precisely because the wave is
* conflict-free: no wave member touches the cells of another, so no member's
* Cell_handle can be destroyed or recycled out from under it while the wave
* runs.
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
  using Edge = typename Operation::Edge;
  using Vertex_handle = typename C3t3::Vertex_handle;
  using Cell_handle = typename C3t3::Cell_handle;
  using FT = typename C3t3::Triangulation::Geom_traits::FT;
  using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
  using Delta = Short_edges_delta<Vertex_pair, FT>;

  // Overridable at run time for the same reason as the other CGAL_TR_* knobs
  // in this codebase (see Elementary_operation.h::buckets_per_thread() /
  // use_unordered_waves()): A/B has to be interleaved within one binary and
  // one thermal state, and that's impossible if the choice needs a rebuild.
  // Read once; default is the parallel refresh (the new path).
  // CGAL_TR_PAR_REFRESH=0 selects the old plain serial loop.
  static bool use_par_refresh()
  {
    static const bool on = []
      {
        const char* const e = std::getenv("CGAL_TR_PAR_REFRESH");
        return e == nullptr || *e != '0';
      }();
    return on;
  }

  // Wave-size knob, overridable at run time (CGAL_TR_COLLAPSE_WAVE_FACTOR).
  // wave_size = factor * nthreads, and the selection scan window is
  // 4 * wave_size -- so the SERIAL selection cost per wave grows linearly
  // with thread count while the number of conflict-free edges a wave can
  // actually hold saturates on the mesh's geometry. That is the suspected
  // mechanism behind wave selection measuring SLOWER at 4 threads than at 1
  // (2.95 s -> 3.23 s, see the "COLLAPSE ANATOMY" log entry). The knob makes
  // the wave size separable from nthreads so the hypothesis can be A/B'd in
  // one binary.
  static std::size_t collapse_wave_factor()
  {
    static const std::size_t f = []() -> std::size_t
      {
        const char* const e = std::getenv("CGAL_TR_COLLAPSE_WAVE_FACTOR");
        if (e != nullptr)
        {
          const long v = std::strtol(e, nullptr, 10);
          if (v > 0) return static_cast<std::size_t>(v);
        }
        return static_cast<std::size_t>(COLLAPSE_WAVE_FACTOR);
      }();
    return f;
  }

public:
  bool execute(Operation& op, C3t3& c3t3) const
  {
#ifdef CGAL_TET_REMESHING_COLLAPSE_SEQUENTIAL
    // Opt out of the conflict-free wave and run collapse sequentially. See the
    // note above the class for when that is the better choice.
    return Elementary_operation_execution_sequential<Operation>().execute(op, c3t3);
#else
    Short_edges short_edges = op.get_elements(c3t3);
    if (short_edges.empty())
      return false;

    ensure_lock_data_structure_initialized(c3t3);

    auto& tr = c3t3.triangulation();
    const std::size_t nthreads =
      static_cast<std::size_t>((std::max)(1, tbb::this_task_arena::max_concurrency()));
    // Collapse wants far larger waves than split's ORDERED_WAVE_FACTOR: only a
    // fraction of the candidates scanned survive the conflict test, so a wave
    // sized like split's leaves most threads idle.
    const std::size_t wave_size =
      (std::max)(std::size_t(1), collapse_wave_factor() * nthreads);

    std::vector<Edge> wave;
    wave.reserve(wave_size);
    std::vector<typename Short_edges::right_map::iterator> selected;
    selected.reserve(wave_size);
    boost::container::small_vector<Cell_handle, 64> inc_cells;
    // surviving vertices whose incident edges the wave deferred to the replay
    std::vector<Vertex_handle> refresh_vertices;
    // drawn from the shared monotonic counter, so stale stamps -- including
    // ones left by the unordered wave executor -- never read as claimed
    std::size_t wave_stamp = 0;
    tbb::combinable<Delta> deltas;
    std::vector<std::pair<Vertex_pair, std::optional<FT>>> replay;
    std::vector<std::size_t> order;
    // per-vertex results of the deferred refresh, evaluated in parallel below;
    // kept outside the wave loop so its inner Delta vectors keep their
    // capacity across waves instead of reallocating every time
    std::vector<Delta> refresh_out;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    CGAL::Real_timer timer, timer_select, timer_par, timer_replay, timer_refresh;
    // Selection-yield counters. `scanned` is the serial work the selection
    // loop does; `accepted` is what it buys. accepted/scanned falling as
    // nthreads rises is the signature of the anti-scaling described above.
    std::size_t n_waves = 0, n_scanned = 0, n_accepted = 0, n_stale = 0;
    timer.start();
#endif

#ifdef CGAL_TETRAHEDRAL_REMESHING_LIVELOCK_GUARD
    // Livelock diagnostic. This loop terminates only by `short_edges` emptying
    // or by the `wave.empty() && selected.empty()` break, and the post-join
    // refresh re-inserts edges around every surviving vertex -- so an edge that
    // is persistently is_too_short() but that collapse_edge() always refuses is
    // erased and re-inserted forever, spinning the main thread while the workers
    // idle. Track the smallest work-list size seen: real progress keeps setting
    // new minima, a livelock stops. On stall, name the offending edge instead of
    // hanging.
    std::size_t best_size = (std::numeric_limits<std::size_t>::max)();
    std::size_t stalled_iters = 0;
    const std::size_t stall_limit = 10000;
#endif
    while (!short_edges.empty())
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_LIVELOCK_GUARD
      if (short_edges.size() < best_size)
      {
        best_size = short_edges.size();
        stalled_iters = 0;
      }
      else if (++stalled_iters > stall_limit)
      {
        const auto it = short_edges.right.begin();
        std::cerr << "[livelock] collapse wave made no progress for "
                  << stall_limit << " iterations.\n"
                  << "  work list size " << short_edges.size()
                  << " (best seen " << best_size << ")\n";
        if (it != short_edges.right.end())
        {
          const Vertex_handle a = it->second.first;
          const Vertex_handle b = it->second.second;
          Cell_handle c; int i0, i1;
          std::cerr << "  head edge sqlen " << it->first
                    << "  dims " << a->in_dimension() << "/" << b->in_dimension()
                    << "  still_an_edge " << tr.tds().is_edge(a, b, c, i0, i1)
                    << "\n  head points " << point(a->point())
                    << " -- " << point(b->point()) << std::endl;
        }
        CGAL_error_msg("tetrahedral remeshing: parallel collapse livelock");
      }
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      timer_select.start();
#endif
      // Select a conflict-free wave: walk the work list shortest-first and take
      // an edge only if no vertex of its incident cells has been claimed by an
      // edge already in the wave. Two selected edges therefore share no cell,
      // and not even a vertex between their cells, so the regions their
      // collapses read and rewrite are disjoint -- the wave needs no locking
      // at all. Rejected edges stay in the work list and are reconsidered by
      // the next wave, so global shortest-first priority is preserved.
      wave.clear();
      selected.clear();
      wave_stamp = next_wave_claim_stamp();
      // Bound the scan. Conflicting edges are left in the work list, and
      // near the front they are mostly neighbours of edges already selected,
      // so scanning further and further down to fill a wave costs far more
      // than it saves -- an unbounded scan re-examines the whole list every
      // wave and dominates the phase. Whatever the window yields is enough;
      // the rest is picked up by the next wave.
      const std::size_t scan_limit = 4 * wave_size;
      std::size_t scanned = 0;
      for (auto eit = short_edges.right.begin();
           eit != short_edges.right.end()
             && wave.size() < wave_size && scanned < scan_limit; ++eit, ++scanned)
      {
        const Vertex_handle v0 = eit->second.first;
        const Vertex_handle v1 = eit->second.second;

        // The conflict test, O(1). It is exact -- not a filter -- because an
        // accepted edge claims its 2-ring (see below), which is precisely the set
        // of vertices whose own region would overlap this one's.
        if (v0->wave_claim_stamp() == wave_stamp
            || v1->wave_claim_stamp() == wave_stamp)
          continue;

        // resolve the stored vertex pair to a current Edge, keeping its
        // orientation; drop the entry if it is no longer an edge
        Cell_handle cell;
        int i0, i1;
        if (!tr.tds().is_edge(v0, v1, cell, i0, i1))
        {
          selected.push_back(eit); // erase it from the work list below
          continue;
        }
        const Edge e(cell, i0, i1);

        // Extra vertices whose full star this candidate's own collapse_edge()
        // will read during the parallel phase, outside v0/v1 -- see
        // topology_test_link_vertices()'s doc comment and the "ROOT CAUSE
        // FOUND" log entry. Empty for interior edges (the common case), so
        // this costs one facet circulation plus O(1) boundary tests per
        // candidate scanned, and only reaches further when that circulation
        // finds a boundary-incident facet.
        const auto link_vs = op.topology_test_link_vertices(e, c3t3);

        bool link_conflict = false;
        for (const Vertex_handle& vi : link_vs)
        {
          if (vi->wave_claim_stamp() == wave_stamp)
          {
            link_conflict = true;
            break;
          }
        }
        if (link_conflict)
          continue; // leave in work list; retry next wave

        inc_cells.clear();
        tr.incident_cells(v0, std::back_inserter(inc_cells));
        tr.incident_cells(v1, std::back_inserter(inc_cells));
        for (const Vertex_handle& vi : link_vs)
          tr.incident_cells(vi, std::back_inserter(inc_cells));

        // Accepted: claim the region -- every vertex of every cell incident to
        // either endpoint, plus (when non-empty) every vertex of every cell
        // incident to a topology_test() link vertex above. Combined with the
        // conflict tests above this makes wave members disjoint not just in
        // the cells collapse_edge() itself mutates, but also in every star its
        // own topology_test() call reads -- which is what the parallel phase
        // actually needs. The endpoint-only claim alone was proven insufficient
        // by the "ROOT CAUSE FOUND" log entry: get_collapse_type()'s
        // topology_test() reads link vertices' stars from inside collapse_edge,
        // during the parallel phase, which the endpoint claim never covered.
        //
        // A full unconditional 2-ring claim was measured earlier as the
        // alternative (66.4-70.7s against 27.7s for the collapse phase) and
        // rejected on cost; the claim above is cheaper because it only widens
        // for the (usually rare) candidates with a boundary-incident facet,
        // instead of every candidate scanned.
        for (const Cell_handle& c : inc_cells)
          for (int k = 0; k < 4; ++k)
            c->vertex(k)->set_wave_claim_stamp(wave_stamp);

        wave.push_back(e);
        selected.push_back(eit);
      }

      // take the selected edges out of the work list
      for (auto& eit : selected)
        short_edges.right.erase(eit);
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      timer_select.stop();
      ++n_waves;
      n_scanned  += scanned;
      n_accepted += wave.size();
      n_stale    += selected.size() - wave.size();
#endif

      if (wave.empty())
      {
        if (selected.empty())
          break; // no progress possible
        continue; // the window held only stale entries, now dropped
      }

      deltas.clear();
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      timer_par.start();
#endif
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, wave.size()),
        [&](const tbb::blocked_range<std::size_t>& r)
        {
          Delta& delta = deltas.local();
          for (std::size_t i = r.begin(); i != r.end(); ++i)
            op.execute_operation(wave[i], c3t3, delta);
        });

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      timer_par.stop();
      timer_replay.start();
#endif

      // Replay the wave's recorded work-list updates into the bimap, but
      // deduplicate first. collapse() removes all six edges of every cell
      // incident to the vertex it deletes, and adjacent cells share edges, so
      // the same key is typically recorded several times per collapse. Each
      // duplicate that survives to the bimap costs a tree descent and a node
      // update, and this replay is serial, so it grew into the phase's
      // bottleneck. Waves are conflict-free, so duplicate keys only ever come
      // from a single collapse; keeping the last write recorded for a key is
      // exactly what applying them all in order would leave behind.
      replay.clear();
      refresh_vertices.clear();
      deltas.combine_each([&](Delta& delta)
        {
          replay.insert(replay.end(), delta.writes.begin(), delta.writes.end());
          refresh_vertices.insert(refresh_vertices.end(),
                                  delta.refresh_vertices.begin(),
                                  delta.refresh_vertices.end());
          delta.clear();
        });

      order.resize(replay.size());
      std::iota(order.begin(), order.end(), std::size_t(0));
      const Compare_vertex_pairs<Vertex_handle> less;
      tbb::parallel_sort(order.begin(), order.end(),
        [&](const std::size_t a, const std::size_t b)
        {
          if (less(replay[a].first, replay[b].first)) return true;
          if (less(replay[b].first, replay[a].first)) return false;
          return a < b; // ties keep the order the writes were recorded in
        });

      for (std::size_t i = 0; i < order.size(); ++i)
      {
        // within a run of equal keys, only the last write survives
        if (i + 1 < order.size()
            && !less(replay[order[i]].first, replay[order[i + 1]].first))
          continue;
        auto& w = replay[order[i]];
        update_bimap(w.first, short_edges, w.second);
      }

      // The wave has joined, so the mesh is stable again : now re-evaluate the
      // edges incident to each surviving vertex. This is the work deferred out
      // of execute_operation (see refresh_incident_edges) because it walks the
      // link vertices' stars, which no wave member claims. Running it here keeps
      // the parallel phase's claim at the cheap endpoint/region level. The
      // ordering matches the sequential path : collapse_edge's own work-list
      // writes are replayed first, then the incident refresh overwrites them.
      //
      // The mesh is stable here, so this re-evaluation is read-only over the
      // TDS : its only side effect is the work-list write. Evaluate every
      // vertex's incident edges in parallel and apply the results serially,
      // in index order, which is exactly the sequence the serial loop
      // produced -- so last-write-wins on duplicate keys is preserved
      // bit-for-bit.
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      timer_refresh.start();
#endif
      if (use_par_refresh())
      {
        refresh_out.resize(refresh_vertices.size());
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, refresh_vertices.size()),
          [&](const tbb::blocked_range<std::size_t>& r)
          {
            for (std::size_t i = r.begin(); i != r.end(); ++i)
            {
              refresh_out[i].clear();
              const Vertex_handle vh = refresh_vertices[i];
              if (!c3t3.triangulation().tds().vertices().is_used(vh))
                continue; // paranoia: a later wave may already have consumed it
              op.refresh_incident_edges(vh, c3t3, refresh_out[i]);
            }
          });
        for (std::size_t i = 0; i < refresh_vertices.size(); ++i)
          for (const auto& w : refresh_out[i].writes)
            update_bimap(w.first, short_edges, w.second);
      }
      else
      {
        // CGAL_TR_PAR_REFRESH=0: old plain serial loop, kept for A/B timing.
        for (const Vertex_handle& vh : refresh_vertices)
        {
          if (!c3t3.triangulation().tds().vertices().is_used(vh))
            continue;
          op.refresh_incident_edges(vh, c3t3, short_edges);
        }
      }
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      timer_refresh.stop();
      timer_replay.stop();
#endif
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    timer.stop();
    std::cout << op.operation_name() << " (parallel): " << timer.time()
              << " sec [wave selection " << timer_select.time()
              << " sec, parallel geometry " << timer_par.time()
              << " sec, serial work-list replay " << timer_replay.time()
              << " sec (of which deferred refresh " << timer_refresh.time() << ")"
              << " sec]." << std::endl;
    std::cout << "  [wave-yield] waves " << n_waves
              << " | wave_size " << wave_size
              << " | scanned " << n_scanned
              << " | accepted " << n_accepted
              << " | stale " << n_stale
              << " | yield " << (n_scanned ? double(n_accepted) / double(n_scanned) : 0.)
              << " | scanned/accepted "
              << (n_accepted ? double(n_scanned) / double(n_accepted) : 0.)
              << std::endl;
#endif

    return true;
#endif
  }

private:
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
};

#endif // CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && CGAL_LINKED_WITH_TBB

}
}
}

#endif // CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H
