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

#ifndef CGAL_INTERNAL_FLIP_EDGES_H
#define CGAL_INTERNAL_FLIP_EDGES_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/utility.h>

#include <CGAL/Tetrahedral_remeshing/internal/Elementary_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>
#include <type_traits>
#include <limits>
#include <queue>
#include <mutex>
#include <algorithm>

#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <cstdlib>
#endif

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{

#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
// Candidate collection for both flip operations used to walk
// Triangulation_3::finite_edges() serially. That iterator pays a
// cell-circulator canonical-owner dedup per edge, and the walk dominates
// get_elements() -- pure serial time inside an otherwise parallel phase, so it
// is bounded by Amdahl rather than by the core count.
// parallel_collect_finite_edges() (tetrahedral_remeshing_helpers.h) applies the
// same canonical-owner rule from a parallel cell scan, and is already what
// split and collapse use for their own candidate collection. Safe here because
// get_elements() runs before the parallel apply pass, with the mesh static, and
// the predicates involved (is_internal / is_boundary) only read.
//
// Runtime toggle rather than compile-time so both arms live in one binary and
// can be interleaved ABBA within a single thermal state.
inline bool par_flip_collect()
{
  static const bool on = []
    {
      const char* const e = std::getenv("CGAL_TR_PAR_FLIP_COLLECT");
      return !(e != nullptr && *e == '0');
    }();
  return on;
}

#ifndef CGAL_TR_FLIP_BUCKETS_PER_THREAD_DEFAULT
// Buckets per thread for the internal flip's kd-partition, overriding the
// global LB_BUCKETS_PER_THREAD (which stays 1, the right value for smooth).
//
// The global default was tuned when "only 15.6% of lock attempts ever retry",
// i.e. when the cost being minimised was the *number of locks acquired* and
// extra buckets only added interface area. That premise no longer holds for
// flip: since the crash fix, the internal flip is the one operation with
// prefers_bucket_scoped_locks() == true, so a worker holds its locks for the
// whole bucket and an interface element cannot acquire until the neighbouring
// bucket completes -- roughly 100k elements of waiting at one bucket/thread.
// Measured on bear at 4 threads, internal flip issues ~12 lock retries per
// locked element while smooth, on the identical executor path, issues 0.03.
// Splitting into many short-lived buckets bounds that wait.
//
// Set via CGAL_TR_FLIP_BUCKETS_PER_THREAD at run time; 0 falls back to the
// global value, which is the A/B baseline arm.
#  define CGAL_TR_FLIP_BUCKETS_PER_THREAD_DEFAULT 32
#endif

inline std::size_t flip_buckets_per_thread()
{
  static const std::size_t n = []() -> std::size_t
    {
      const char* const e = std::getenv("CGAL_TR_FLIP_BUCKETS_PER_THREAD");
      if (e != nullptr)
      {
        const long v = std::strtol(e, nullptr, 10);
        if (v >= 0)
          return static_cast<std::size_t>(v);
      }
      return static_cast<std::size_t>(CGAL_TR_FLIP_BUCKETS_PER_THREAD_DEFAULT);
    }();
  return n;
}
#endif

enum Flip_Criterion{ MIN_ANGLE_BASED, AVERAGE_ANGLE_BASED,
                     VALENCE_BASED, VALENCE_MIN_DH_BASED };

// ---------------------------------------------------------------------------
// [DEBUG-lk01] Lock-ownership probe.
//
// Diagnostic only: enable with -DCGAL_TR_FLIP_LOCK_PROBE.
//
// The internal-flip crash is a read of the cell star of a ring vertex `w` that
// races a concurrent flip mutating those same cells. That race is timing
// dependent and neither gdb nor TSan reproduce it (both serialise the run hard
// enough that the workers never overlap). This probe sidesteps timing entirely:
// instead of waiting to observe two accesses collide, it asserts the invariant
// that would make a collision impossible -- "before I walk w's cell star, I
// hold w's spatial lock". A violation is reported on the very first run that
// executes the offending code path, whether or not another worker happened to
// be in the window at that moment.
//
// Reports (does not abort) so a single run enumerates every distinct offending
// site rather than stopping at the first.
// ---------------------------------------------------------------------------
#ifdef CGAL_TR_FLIP_LOCK_PROBE
template<typename Tr, typename Vertex_handle>
void tr_flip_lock_probe(const Tr& tr, Vertex_handle v, const char* site)
{
  using Lds = std::remove_pointer_t<decltype(tr.get_lock_data_structure())>;
  if constexpr (std::is_void_v<Lds>)
    return; // sequential triangulation: no lock grid to check against
  else
  {
  auto* lds = tr.get_lock_data_structure();
  if (lds == nullptr)
    return; // locking disabled: nothing to check
  if (lds->is_locked_by_this_thread(v->point()))
    return;

  const int mode = tr_probe_mode();
  const char* mode_name = (mode == 1) ? "INTERIOR-elided"
                        : (mode == 2) ? "LOCKED-ZONE"
                        : (mode == 3) ? "WAVE"
                                      : "other/serial";

  // Rate-limit per (site, mode) so one hot loop cannot flood the log, but a
  // rare LOCKED-ZONE violation is never hidden behind a common INTERIOR one.
  static std::mutex probe_mutex;
  static std::unordered_map<std::string, std::size_t> hits;
  std::lock_guard<std::mutex> guard(probe_mutex);
  const std::size_t n = ++hits[std::string(mode_name) + "|" + site];
  if (n <= 3 || n % 5000 == 0)
  {
    std::cerr << "[DEBUG-lk01] UNLOCKED READ  mode=" << mode_name
              << "  site=" << site
              << "  (hit #" << n << ")"
              << "  vertex=" << (&*v)
              << "  locked_by_someone=" << lds->is_locked(v->point())
              << std::endl;
  }
  }
}
#define CGAL_TR_LOCK_PROBE(tr, v, site) \
  ::CGAL::Tetrahedral_remeshing::internal::tr_flip_lock_probe((tr), (v), (site))
#else
#define CGAL_TR_LOCK_PROBE(tr, v, site) CGAL_USE(tr)
#endif

// ---------------------------------------------------------------------------
// Lock-on-demand for the flip footprint.
//
// The precomputed-footprint approach (Elementary_operation::locked_vertices(),
// evaluated in Pass 1) cannot protect a flip: it is computed before any
// parallel mutation has happened, but a flip's actual footprint at execute
// time depends on the neighbourhood as it is *then*, which concurrent flips
// have already reshaped. Widening that precomputed set does not help, because
// the set is computed at the wrong time rather than merely being too small.
//
// Instead, every read of a vertex's cell star acquires that vertex's lock at
// the moment of the read. If the lock is unavailable, the whole flip is
// abandoned via Flip_lock_bail: flips are best-effort (the edge simply stays
// unflipped and may be revisited on a later pass), so abandoning one costs
// nothing but a little work. Crucially, all such reads happen during flip
// *selection*, before any cell is created or deleted, so a bail can never
// leave the triangulation partially modified.
//
// try_lock_vertex() never blocks, and the executor calls
// unlock_all_elements() after each element, so this cannot deadlock.
// ---------------------------------------------------------------------------
struct Flip_lock_bail {};

template<typename Tr, typename Vertex_handle>
inline void tr_flip_require_lock(const Tr& tr, Vertex_handle v)
{
#ifdef CGAL_TR_FLIP_NO_LOCK_ON_DEMAND
  // A/B switch: restores the pre-fix behaviour (precomputed footprint only)
  // so the cost of lock-on-demand can be measured against it. Crashes.
  CGAL_USE(tr); CGAL_USE(v);
  return;
#else
  // A sequential triangulation has no lock grid at all and declares
  // get_lock_data_structure() as returning void*, so the body below cannot
  // even be compiled for it. Nothing can race us there either.
  using Lds = std::remove_pointer_t<decltype(tr.get_lock_data_structure())>;
  if constexpr (std::is_void_v<Lds>)
  {
    CGAL_USE(tr); CGAL_USE(v);
  }
  else
  {
    auto* lds = tr.get_lock_data_structure();
    if (lds == nullptr)
      return; // parallel build, but locking disabled for this run

    if (lds->is_locked_by_this_thread(v->point()))
      return;

    if (!const_cast<Tr&>(tr).try_lock_vertex(v))
      throw Flip_lock_bail{};
  }
#endif
}
#define CGAL_TR_REQUIRE_LOCK(tr, v) \
  ::CGAL::Tetrahedral_remeshing::internal::tr_flip_require_lock((tr), (v))

//outer_mirror_facets contains the set of facets of the outer hull
//of the set of cells modified by the flip operation,
//"seen from" outside
//i.e. for each facet f among those, f.first has not been modified by flip
template<typename C3t3, typename CellSet, typename FacetSet>
void update_c3t3_facets(C3t3& c3t3,
                        const CellSet& cells_to_update,
                        const FacetSet& outer_mirror_facets)
{
  typedef typename C3t3::Facet       Facet;
  typedef typename C3t3::Cell_handle Cell_handle;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;

  for (Cell_handle c : cells_to_update)
  {
    //their subdomain indices have not been modified because we kept the same cells
    //surface patch indices need to be fixed though
    for (int i = 0; i < 4; ++i)
    {
      const Facet f(c, i);
      const Facet mf = c3t3.triangulation().mirror_facet(f);
      if (outer_mirror_facets.find(mf) != outer_mirror_facets.end())
      {
        //we are on the border of the modified zone, c3t3 info is valid outside,
        //on mirror facet
        const typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(mf);
        if (c3t3.is_in_complex(mf))
          f.first->set_surface_patch_index(f.second, patch);
        else
          f.first->set_surface_patch_index(f.second, Surface_patch_index());
      }
      else
      {
        //we are inside the modified zone, c3t3 info is not valid anymore
        if (c3t3.is_in_complex(f) || c3t3.is_in_complex(mf))
        {
          f.first->set_surface_patch_index(f.second, Surface_patch_index());
          mf.first->set_surface_patch_index(mf.second, Surface_patch_index());
        }
      }
    }
  }
}

template<typename C3t3, typename IncCellsVectorMap, typename CellSelector>
Sliver_removal_result flip_3_to_2(typename C3t3::Edge& edge,
                                  C3t3& c3t3,
                                  const std::vector<typename C3t3::Vertex_handle>& vertices_around_edge,
                                  const Flip_Criterion& criterion,
                                  IncCellsVectorMap& inc_cells,
                                  CellSelector& cell_selector)
{
  typedef typename C3t3::Triangulation Tr;
  typedef typename C3t3::Facet         Facet;
  typedef typename C3t3::Vertex_handle Vertex_handle;
  typedef typename C3t3::Cell_handle   Cell_handle;
  typedef typename Tr::Cell_circulator Cell_circulator;
  typedef typename Tr::Geom_traits     Gt;
  typedef typename Gt::FT              FT;

  //Edge to face flip
  Tr& tr = c3t3.triangulation();

  Cell_circulator circ = tr.incident_cells(edge);
  Cell_circulator done = circ;

  Vertex_handle vh0 = edge.first->vertex(edge.second);
  Vertex_handle vh1 = edge.first->vertex(edge.third);

  //Select 2 cells to keep and update and one to remove
  Cell_handle ch0 = Cell_handle(circ++);
  Cell_handle ch1 = Cell_handle(circ++);
  Cell_handle cell_to_remove = Cell_handle(circ++);
  if (circ != done)
  {
    std::cout << "Wrong flip function" << std::endl;
    return NOT_FLIPPABLE;
  }

  //Check structural validity
  Cell_handle c;
  int i0, i1, i3;
  if (tr.is_facet(vertices_around_edge[0], vertices_around_edge[1], vertices_around_edge[2],
                  c, i0, i1, i3))
    return NOT_FLIPPABLE;

  //Check topological validity
  const typename C3t3::Subdomain_index subdomain = ch0->subdomain_index();
  if ( subdomain != ch1->subdomain_index()
       || subdomain != cell_to_remove->subdomain_index()
       || ch1->subdomain_index() != cell_to_remove->subdomain_index())
    return NOT_FLIPPABLE;

  Vertex_handle vh2;
  Vertex_handle vh3;

  for (int i = 0; i < 3; ++i){
    if (!ch0->has_vertex(vertices_around_edge[i]))
      vh2 = vertices_around_edge[i];
    else if (!ch1->has_vertex(vertices_around_edge[i]))
      vh3 = vertices_around_edge[i];
  }

  int vh0_id = ch0->index(vh0);
  int vh1_id = ch1->index(vh1);

  //Check if flip valid
  if (!is_well_oriented(tr, vh2,
                        ch0->vertex(indices(vh0_id, 0)),
                        ch0->vertex(indices(vh0_id, 1)),
                        ch0->vertex(indices(vh0_id, 2)))
      || !is_well_oriented(tr, vh3,
                           ch1->vertex(indices(vh1_id, 0)),
                           ch1->vertex(indices(vh1_id, 1)),
                           ch1->vertex(indices(vh1_id, 2))))
    return NOT_FLIPPABLE;

  ///********************VALIDITY CHECK***************************/
  //double curr_min_dh;
  //bool check_validity = false;
  //std::vector<typename Tr::Tetrahedron> pre_sliver_Removal_cells;
  //if (check_validity){
  //  pre_sliver_Removal_cells.clear();
  //  pre_sliver_Removal_cells.push_back(K::Tetrahedron_3(ch0->vertex(0)->point(), ch0->vertex(1)->point(), ch0->vertex(2)->point(), ch0->vertex(3)->point()));
  //  pre_sliver_Removal_cells.push_back(K::Tetrahedron_3(ch1->vertex(0)->point(), ch1->vertex(1)->point(), ch1->vertex(2)->point(), ch1->vertex(3)->point()));
  //  pre_sliver_Removal_cells.push_back(K::Tetrahedron_3(cell_to_remove->vertex(0)->point(), cell_to_remove->vertex(1)->point(),
  //    cell_to_remove->vertex(2)->point(), cell_to_remove->vertex(3)->point()));

  //  curr_min_dh = min_dihedral_angle<Gt>(ch0);
  //  curr_min_dh = std::min(curr_min_dh, min_dihedral_angle<Gt>(ch1));
  //  curr_min_dh = std::min(curr_min_dh, min_dihedral_angle<Gt>(cell_to_remove));

  //  pre_sliver_Removal_vertices.clear();
  //  for (int i = 0; i < vertices_around_edge.size(); ++i){
  //    pre_sliver_Removal_vertices.push_back(Point_3(vertices_around_edge[i]->point()));
  //  }

  //  previous_edges.clear();
  //  previous_edges.push_back(std::make_pair(vh0->point(), vh1->point()));
  //}
  /*************************************************************/


  if (criterion == MIN_ANGLE_BASED)
  {
    //Current worst dihedral angle
    Dihedral_angle_cosine curr_max_cosdh = max_cos_dihedral_angle(tr, ch0);
    curr_max_cosdh = (std::max)(curr_max_cosdh, max_cos_dihedral_angle(tr, ch1));
    curr_max_cosdh = (std::max)(curr_max_cosdh, max_cos_dihedral_angle(tr, cell_to_remove));

    //Result worst dihedral angle
    if (curr_max_cosdh < max_cos_dihedral_angle(tr, vh2,
                                         ch0->vertex(indices(vh0_id, 0)),
                                         ch0->vertex(indices(vh0_id, 1)),
                                         ch0->vertex(indices(vh0_id, 2)))
        || curr_max_cosdh < max_cos_dihedral_angle(tr, vh3,
                                         ch1->vertex(indices(vh1_id, 0)),
                                         ch1->vertex(indices(vh1_id, 1)),
                                         ch1->vertex(indices(vh1_id, 2))))
      return NO_BEST_CONFIGURATION;
  }
  else if (criterion == AVERAGE_ANGLE_BASED)
  {
    //Current worst dihedral angle
    double average_min_dh = min_dihedral_angle(tr, ch0);
    average_min_dh += min_dihedral_angle(tr, ch1);
    average_min_dh += min_dihedral_angle(tr, cell_to_remove);

    average_min_dh /= 3.;

    FT new_average_min_dh = 0.5 *
                            (min_dihedral_angle(tr, vh2, ch0->vertex(indices(vh0_id, 0)),
                                ch0->vertex(indices(vh0_id, 1)),
                                ch0->vertex(indices(vh0_id, 2)))
                           + min_dihedral_angle(tr, vh3, ch1->vertex(indices(vh1_id, 0)),
                                 ch1->vertex(indices(vh1_id, 1)),
                                 ch1->vertex(indices(vh1_id, 2))));
    //Result worst dihedral angle
    if (average_min_dh > new_average_min_dh)
      return NO_BEST_CONFIGURATION;
  }

  //Keep the facets
  typedef CGAL::Triple<Vertex_handle, Vertex_handle, Vertex_handle> Facet_vvv;
  typedef std::unordered_map<Facet_vvv, std::size_t> FaceMapIndex;
  std::unordered_set<Facet, boost::hash<Facet>> outer_mirror_facets;

  FaceMapIndex facet_map_indices;
  std::vector<Facet> mirror_facets;
  circ = Cell_circulator(done);
  do
  {
    // facet opposite to vh0
    int curr_vh0_id = circ->index(vh0);
    Facet n_vh0_facet = tr.mirror_facet(Facet(circ, curr_vh0_id));

    outer_mirror_facets.insert(n_vh0_facet);

    Facet_vvv face0 = make_vertex_triple(circ->vertex(indices(curr_vh0_id, 0)),
                                         circ->vertex(indices(curr_vh0_id, 1)),
                                         circ->vertex(indices(curr_vh0_id, 2)));

    typename FaceMapIndex::iterator it = facet_map_indices.find(face0);
    if (it == facet_map_indices.end())
    {
      facet_map_indices[face0] = mirror_facets.size();
      mirror_facets.push_back(n_vh0_facet);
    }

    // facet opposite to vh1
    int curr_vh1_id = circ->index(vh1);
    Facet n_vh1_facet = tr.mirror_facet(Facet(circ, curr_vh1_id));

    outer_mirror_facets.insert(n_vh1_facet);

    Facet_vvv face1 = make_vertex_triple(circ->vertex(indices(curr_vh1_id, 0)),
                                         circ->vertex(indices(curr_vh1_id, 1)),
                                         circ->vertex(indices(curr_vh1_id, 2)));
    it = facet_map_indices.find(face1);
    if (it == facet_map_indices.end())
    {
      facet_map_indices[face1] = mirror_facets.size();
      mirror_facets.push_back(n_vh1_facet);
    }
  }
  while (++circ != done);

  /*
  c3t3.remove_from_complex( ch0 );
  c3t3.remove_from_complex( ch1 );
  c3t3.remove_from_complex( cell_to_remove );

  tr.flip(edge);

  for( int i = 0 ; i < facets.size() ; i ++ ){
  Cell_handle new_cell = facets[i].first->neighbor( facets[i].second );
  c3t3.add_to_complex( new_cell, si );
  }
  */

  //Update cells
  ch0->set_vertex(vh0_id, vh2);
  ch1->set_vertex(vh1_id, vh3);

  // "New" cells are not created, only modified/updated
  std::vector<Cell_handle> cells_to_update;
  cells_to_update.push_back(ch0);
  cells_to_update.push_back(ch1);

  //Update adjacencies and vertices' cells
  for (Cell_handle ch : cells_to_update)
  {
    for (int v = 0; v < 4; ++v)
    {
      Facet_vvv face = make_vertex_triple(ch->vertex(indices(v, 0)),
                                          ch->vertex(indices(v, 1)),
                                          ch->vertex(indices(v, 2)));
      typename FaceMapIndex::iterator it = facet_map_indices.find(face);
      if (it == facet_map_indices.end())
      {
        facet_map_indices[face] = mirror_facets.size();
        mirror_facets.push_back(Facet(ch, v));
      }
      else
      {
        Facet mirror_facet = mirror_facets[it->second];

        //Update neighbor
        mirror_facet.first->set_neighbor(mirror_facet.second, ch);
        ch->set_neighbor(v, mirror_facet.first);
      }
      ch->vertex(v)->set_cell(ch);

      inc_cells[ch->vertex(v)].clear();
      ch->reset_cache_validity();
    }
  }

  // Update c3t3
  update_c3t3_facets(c3t3, cells_to_update, outer_mirror_facets);

  treat_before_delete(cell_to_remove, cell_selector, c3t3);
  tr.tds().delete_cell(cell_to_remove);

  /********************VALIDITY CHECK***************************/
  //if (check_validity)
  //{
  //  post_sliver_Removal_cells.clear();
  //  post_sliver_Removal_cells.push_back(ch0);
  //  post_sliver_Removal_cells.push_back(ch1);

  //  double new_min_dh = min_dihedral_angle<Gt>(ch0);
  //  new_min_dh = std::min(new_min_dh, min_dihedral_angle<Gt>(ch1));

  //  post_sliver_Removal_vertices.clear();
  //  post_sliver_Removal_vertices.push_back(vh2);
  //  post_sliver_Removal_vertices.push_back(vh3);

  //  if (!is_well_oriented(ch0))
  //    return INVALID_ORIENTATION;
  //  if (!is_well_oriented(ch1))
  //    return INVALID_ORIENTATION;
  //  if (!tr.is_valid(ch0))
  //    return INVALID_CELL;
  //  if (!tr.is_valid(ch1))
  //    return INVALID_CELL;

  //  for (int i = 0; i < 4; ++i){
  //    if (!tr.is_valid(ch0->neighbor(i)))
  //      return INVALID_CELL;
  //    if (!tr.is_valid(ch1->neighbor(i)))
  //      return INVALID_CELL;
  //    if (!tr.tds().is_valid(ch0->vertex(i)))
  //      return INVALID_VERTEX;
  //    if (!tr.tds().is_valid(ch1->vertex(i))){
  //      return INVALID_VERTEX;
  //    }
  //  }

  //  if ((curr_min_dh - new_min_dh) > 0.01){
  //    std::cout << "Three_to_two_flip::Flip not improving the quality: " << curr_min_dh << " to " << new_min_dh << std::endl;
  //    return INVALID_CELL;
  //  }
  //}
  /***********************************************************/

  return VALID_FLIP;
}

template<typename C3t3, typename CandidatesQueue>
void find_best_flip_to_improve_dh(C3t3& c3t3,
                                  typename C3t3::Edge& edge,
                                  typename C3t3::Vertex_handle vh2,
                                  typename C3t3::Vertex_handle vh3,
                                  CandidatesQueue& candidates,
                                  const Dihedral_angle_cosine& curr_max_cos_dh,
                                  bool is_sliver_well_oriented = true,
                                  int e_id = 0)
{
  typedef typename C3t3::Triangulation  Tr;
  typedef typename C3t3::Vertex_handle  Vertex_handle;
  typedef typename C3t3::Facet          Facet;
  typedef typename Tr::Facet_circulator Facet_circulator;
  typedef typename Tr::Cell_circulator  Cell_circulator;

  // std::cout << "find_best_flip_to_improve_dh boundary " << std::endl;
  Tr& tr = c3t3.triangulation();

  Vertex_handle vh0 = edge.first->vertex(edge.second);
  Vertex_handle vh1 = edge.first->vertex(edge.third);

  Facet_circulator curr_fcirc = tr.incident_facets(edge);
  Facet_circulator curr_fdone = curr_fcirc;

  //Only keep the possible flips
  std::vector<Vertex_handle> opposite_vertices;
  int nb_cells_around_edge = 0;
  do
  {
    Vertex_handle vh;
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle curr_vertex = curr_fcirc->first->vertex(indices(curr_fcirc->second, i));
      if ( curr_vertex != vh0
           && curr_vertex != vh1
           && (curr_vertex == vh2 || curr_vertex == vh3))
      {
        vh = curr_vertex;
        Facet_circulator facet_circulator(curr_fcirc);
        Facet_circulator facet_done(curr_fcirc);

        facet_done--;
        facet_circulator++;
        facet_circulator++;

        bool is_edge = false;
        do
        {
          //Get the ids of the opposite vertices
          for (int j = 0; j < 3; ++j)
          {
            Vertex_handle curr = facet_circulator->first->vertex(
                                          indices(facet_circulator->second, j));
            if (curr != vh0  && curr != vh1)
            {
              if (tr.tds().is_edge(curr, vh))
                is_edge = true;
            }
          }
        } while (++facet_circulator != facet_done);

        if (!is_edge && !tr.is_infinite(vh))
          opposite_vertices.push_back(vh);
      }
    }
    nb_cells_around_edge++;
  }
  while (++curr_fcirc != curr_fdone);

  if (nb_cells_around_edge < 4)
    return;

  //Facets that will be used to create new cells i.e. all the facets opposite to vh1 and don't have vh
  //Facets that will be used to update cells i.e. all the facets opposite to vh0 will be set to vh: facet.first->set_vertex( facet.second, vh )

  Cell_circulator cell_circulator = tr.incident_cells(edge);
  Cell_circulator done = cell_circulator;

  for (std::size_t i = 0; i < opposite_vertices.size(); ++i)
  {
    Vertex_handle vh = opposite_vertices[i];
    bool keep = true;

    boost::container::small_vector<Facet, 60> facets;
    do
    {
      //Store it if it do not have vh
      if (!cell_circulator->has_vertex(vh))
      {
        //Facets opposite to vh0
        Facet facet_vh0(cell_circulator, cell_circulator->index(vh0));

        //Facets opposite to vh1
        Facet facet_vh1(cell_circulator, cell_circulator->index(vh1));

        facets.push_back(facet_vh1);
        facets.push_back(facet_vh0);
      }
    } while (++cell_circulator != done);


    Dihedral_angle_cosine max_flip_cos_dh(CGAL::NEGATIVE, 1., 1.);
    for (const Facet& fi : facets)
    {
      if (!tr.is_infinite(fi.first) && c3t3.is_in_complex(fi.first))
      {
        if (is_well_oriented(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                             fi.first->vertex(indices(fi.second, 1)),
                             fi.first->vertex(indices(fi.second, 2))))
        {
          max_flip_cos_dh = (std::max)(
            max_flip_cos_dh,
            max_cos_dihedral_angle(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                                           fi.first->vertex(indices(fi.second, 1)),
                                           fi.first->vertex(indices(fi.second, 2))));
        }
        else
        {
          keep = false;
          break;
        }

        if (max_flip_cos_dh.is_one())//it will not get worse than 1.
        {
          keep = false;
          break;
        }

        // the worst angle of the flip only ever grows from here, so once it has
        // reached the one the edge already has, this vertex cannot be kept -
        // unless the cells are inverted, where it is kept whatever it measures
        if (is_sliver_well_oriented && !(max_flip_cos_dh < curr_max_cos_dh))
        {
          keep = false;
          break;
        }
      }
    }
    facets.clear();

    if (keep && (max_flip_cos_dh < curr_max_cos_dh  || !is_sliver_well_oriented))
    {
      //std::cout << "vh " << vh->info() <<" old " << curr_max_cos_dh << " min " << min_flip_tan_dh << std::endl;
      candidates.push(std::make_pair(max_flip_cos_dh, std::make_pair(vh, e_id)));
    }
  }
}

template<typename Vertex_handle, typename CellVector, typename Cell_handle>
bool is_edge_uv(Vertex_handle u,
                Vertex_handle v,
                const CellVector& cells_incident_to_u,
                Cell_handle& cell,
                int& i,
                int& j)
{
  if (u == v)
    return false;

  for (typename CellVector::value_type c : cells_incident_to_u)
  {
    if (c->has_vertex(v, j))
    {
      cell = c;
      i = cell->index(u);
      return true;
    }
  }
  return false;
}

template<typename Vertex_handle, typename CellVector>
bool is_edge_uv(Vertex_handle u,
                Vertex_handle v,
                const CellVector& cells_incident_to_u)
{
  typename CellVector::value_type c;
  int i, j;
  return is_edge_uv(u, v, cells_incident_to_u, c, i, j);
}

template<typename C3t3, typename CandidatesQueue,
         typename IncCellsVectorMap>
void find_best_flip_to_improve_dh(C3t3& c3t3,
                                  typename C3t3::Edge& edge,
                                  CandidatesQueue& candidates,
                                  const Dihedral_angle_cosine& curr_max_cosdh,
                                  IncCellsVectorMap& inc_cells,
                                  bool is_sliver_well_oriented = true,
                                  int e_id = 0)
{
  typedef typename C3t3::Triangulation  Tr;
  typedef typename C3t3::Vertex_handle  Vertex_handle;
  typedef typename C3t3::Cell_handle    Cell_handle;
  typedef typename C3t3::Facet          Facet;
  typedef typename Tr::Facet_circulator Facet_circulator;
  typedef typename Tr::Cell_circulator  Cell_circulator;

  Tr& tr = c3t3.triangulation();

  Vertex_handle vh0 = edge.first->vertex(edge.second);
  Vertex_handle vh1 = edge.first->vertex(edge.third);

  Facet_circulator curr_fcirc = tr.incident_facets(edge);
  Facet_circulator curr_fdone = curr_fcirc;

  //Only keep the possible flips
  std::vector<Vertex_handle> opposite_vertices;
  int nb_cells_around_edge = 0;
  do
  {
    Vertex_handle vh;
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle curr_vertex = curr_fcirc->first->vertex(
                                    indices(curr_fcirc->second, i));
      if (curr_vertex != vh0 && curr_vertex != vh1)
      {
        vh = curr_vertex;
        break;
      }
    }

    if(tr.is_infinite(vh))
      continue;

    boost::container::small_vector<Cell_handle, 64>& o_inc_vh = inc_cells[vh];
    if (o_inc_vh.empty())
    {
      CGAL_TR_LOCK_PROBE(tr, vh, "find_best_flip_to_improve_dh:ring_vh");
      CGAL_TR_REQUIRE_LOCK(tr, vh);
      tr.incident_cells(vh, std::back_inserter(o_inc_vh));
    }

    Facet_circulator facet_circulator = curr_fcirc;
    Facet_circulator facet_done = curr_fcirc;

    facet_done--;
    facet_circulator++;
    facet_circulator++;
    bool is_edge = false;
    do
    {
      //Get the ids of the opposite vertices
      for (int i = 0; i < 3; ++i)
      {
        Vertex_handle curr_vertex = facet_circulator->first->vertex(
                                      indices(facet_circulator->second, i));
        if (curr_vertex != vh0  && curr_vertex != vh1)
        {
          if (is_edge_uv(vh, curr_vertex, o_inc_vh))
          {
            is_edge = true;
            break;
          }
        }
      }
    } while (++facet_circulator != facet_done);

    if (!is_edge)
      opposite_vertices.push_back(vh);

    nb_cells_around_edge++;
  }
  while (++curr_fcirc != curr_fdone);
  if (nb_cells_around_edge < 4)
    return;

  //Facets that will be used to create new cells
  //    i.e. all the facets opposite to vh1 and don't have vh
  //Facets that will be used to update cells
  //    i.e. all the facets opposite to vh0 will be set to vh:
  //    facet.first->set_vertex( facet.second, vh )

  Cell_circulator cell_circulator = tr.incident_cells(edge);
  Cell_circulator done = cell_circulator;

  boost::container::small_vector<Facet, 60> facets;
  for (Vertex_handle vh : opposite_vertices)
  {
    bool keep = true;
    do
    {
      //Store it if it do not have vh
      if (!cell_circulator->has_vertex(vh))
      {
        //Facets opposite to vh0
        Facet facet_vh0(cell_circulator, cell_circulator->index(vh0));

        //Facets opposite to vh1
        Facet facet_vh1(cell_circulator, cell_circulator->index(vh1));

        facets.push_back(facet_vh1);
        facets.push_back(facet_vh0);
      }
    }
    while (++cell_circulator != done);

    Dihedral_angle_cosine max_flip_cos_dh(CGAL::NEGATIVE, 1., 1.);
    for (const Facet& fi : facets)
    {
      if (!tr.is_infinite(fi.first))
      {
        if (is_well_oriented(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                             fi.first->vertex(indices(fi.second, 1)),
                             fi.first->vertex(indices(fi.second, 2))))
        {
          max_flip_cos_dh = (std::max)(max_flip_cos_dh,
            max_cos_dihedral_angle(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                                           fi.first->vertex(indices(fi.second, 1)),
                                           fi.first->vertex(indices(fi.second, 2))));
        }
        else
        {
          keep = false;
          break;
        }

        if (max_flip_cos_dh.is_one())//it will not get worse than 1.
        {
          keep = false;
          break;
        }

        // the worst angle of the flip only ever grows from here, so once it has
        // reached the one the edge already has, this vertex cannot be kept -
        // unless the cells are inverted, where it is kept whatever it measures
        if (is_sliver_well_oriented && !(max_flip_cos_dh < curr_max_cosdh))
        {
          keep = false;
          break;
        }
      }
    }
    facets.clear();

    if (keep && (max_flip_cos_dh < curr_max_cosdh || !is_sliver_well_oriented))
    {
      //std::cout << "vh " << vh->info() <<" old " << curr_max_cosdh << " min " << min_flip_tan_dh << std::endl;
      candidates.push(std::make_pair(max_flip_cos_dh, std::make_pair(vh, e_id)));
    }
  }
}

template<typename C3t3,
         typename IncCellsVectorMap,
         typename CellSelector,
         typename Visitor>
Sliver_removal_result flip_n_to_m(C3t3& c3t3,
                                  typename C3t3::Edge& edge,
                                  typename C3t3::Vertex_handle vh,
                                  IncCellsVectorMap& inc_cells,
                                  CellSelector& cell_selector,
                                  Visitor& visitor,
                                  bool check_validity = false)
{
  CGAL_USE(check_validity);
  // std::cout << "n_to_m_flip::start" << std::endl;
  typedef typename C3t3::Triangulation  Tr;
  typedef typename C3t3::Vertex_handle  Vertex_handle;
  typedef typename C3t3::Cell_handle    Cell_handle;
  typedef typename C3t3::Facet          Facet;
  typedef typename Tr::Facet_circulator Facet_circulator;
  typedef typename Tr::Cell_circulator  Cell_circulator;

  Tr& tr = c3t3.triangulation();

  const Vertex_handle vh0 = edge.first->vertex(edge.second);
  const Vertex_handle vh1 = edge.first->vertex(edge.third);

  //This vertex will have its valence augmenting a lot,
  //TODO take the best one

  //TODO!!!! Check that the created edges do not exist!!!

  boost::container::small_vector<Facet, 2> facets_in_complex;

  Facet_circulator facet_circulator = tr.incident_facets(edge);
  Facet_circulator done_facet_circulator = facet_circulator;
  bool look_for_vh_iterator = true;
  do
  {
    if (c3t3.is_in_complex(*facet_circulator))
    {
      facets_in_complex.push_back(*facet_circulator);
    }

    facet_circulator++;

    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      if (facet_circulator->first->vertex(indices(facet_circulator->second, i)) == vh)
        look_for_vh_iterator = false;
    }

  } while (facet_circulator != done_facet_circulator && look_for_vh_iterator);

  if (look_for_vh_iterator) {
    std::cout << "Vertex not an opposite of the edge!!" << std::endl;
    return NOT_FLIPPABLE;
  }

  Facet_circulator facet_done(facet_circulator);
  facet_done--;
  facet_circulator++;
  facet_circulator++;

  boost::container::small_vector<Cell_handle, 64>& o_inc_vh = inc_cells[vh];
  if (o_inc_vh.empty())
  {
    CGAL_TR_LOCK_PROBE(tr, vh, "flip_n_to_m:vh");
    CGAL_TR_REQUIRE_LOCK(tr, vh);
    tr.incident_cells(vh, std::back_inserter(o_inc_vh));
  }

  do
  {
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle curr_vertex = facet_circulator->first->vertex(
                                    indices(facet_circulator->second, i));
      if (curr_vertex != vh0  && curr_vertex != vh1)
      {
        if (is_edge_uv(vh, curr_vertex, o_inc_vh))
          return NOT_FLIPPABLE;
      }
    }
  } while (++facet_circulator != facet_done);


  boost::container::small_vector<Cell_handle, 20> to_remove;

  //Neighbors that will need to be updated after flip
  std::unordered_set<Facet, boost::hash<Facet>> neighbor_facets;

  //Facets that will be used to create new cells
  // i.e. all the facets opposite to vh1 and don't have vh
  std::vector<Facet> facets_for_new_cells;

  //Facets that will be used to update cells
  // i.e. all the facets opposite to vh0 will be set to vh :
  // facet.first->set_vertex( facet.second, vh )
  std::vector<Facet> facets_for_updated_cells;

  Cell_circulator cell_circulator = tr.incident_cells(edge);
  Cell_circulator done = cell_circulator;
  do
  {
    //Facets opposite to vh0
    Facet facet_vh0(cell_circulator, cell_circulator->index(vh0));
    neighbor_facets.insert(tr.mirror_facet(facet_vh0));

    //Facets opposite to vh1
    Facet facet_vh1(cell_circulator, cell_circulator->index(vh1));
    neighbor_facets.insert(tr.mirror_facet(facet_vh1));

    //Store it if it do not have vh
    if (cell_circulator->has_vertex(vh)) {
      to_remove.push_back(cell_circulator);
    }
    else
    {
      facets_for_new_cells.push_back(facet_vh1);
      facets_for_updated_cells.push_back(facet_vh0);
    }
    //
    //        if( ! is_well_oriented( cell_circulator ) )
    //            return WRONG;
  }
  while (++cell_circulator != done);

  //Check that the result will be valid
  for (const Facet& fi : facets_for_new_cells)
  {
    if ( !tr.is_infinite(fi.first)
         && !is_well_oriented(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                              fi.first->vertex(indices(fi.second, 1)),
                              fi.first->vertex(indices(fi.second, 2))))
      return NOT_FLIPPABLE;
  }
  for (const Facet& fi : facets_for_updated_cells)
  {
    if ( !tr.is_infinite(fi.first)
         && !is_well_oriented(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                              fi.first->vertex(indices(fi.second, 1)),
                              fi.first->vertex(indices(fi.second, 2))))
      return NOT_FLIPPABLE;
  }

  ///********************VALIDITY CHECK***************************/
  //double current_min_dh = DBL_MAX;

  //if (check_validity){

  //  pre_sliver_Removal_cells.clear();
  //  do{
  //    pre_sliver_Removal_cells.push_back(K::Tetrahedron_3(cell_circulator->vertex(0)->point(), cell_circulator->vertex(1)->point(),
  //      cell_circulator->vertex(2)->point(), cell_circulator->vertex(3)->point()));

  //    if (!tr.is_infinite(cell_circulator))
  //      current_min_dh = std::min(current_min_dh, min_dihedral_angle(cell_circulator));
  //  } while (++cell_circulator != done);

  //  pre_sliver_Removal_vertices.clear();
  //  pre_sliver_Removal_vertices.push_back(vh->point());

  //  previous_edges.clear();
  //  previous_edges.push_back(std::make_pair(vh0->point(), vh1->point()));
  //}
  ///*************************************************************/

  //Surface
  for (const Facet& f : facets_in_complex)
    c3t3.remove_from_complex(f);

  //Subdomain index
  typedef typename C3t3::Subdomain_index Subdomain_index;
  const Subdomain_index subdomain = to_remove[0]->subdomain_index();
  bool selected = get(cell_selector, to_remove[0]);
  visitor.before_flip(to_remove[0]);

  std::vector<Cell_handle> cells_to_update;

  //Create new cells
  for (const Facet& fi : facets_for_new_cells)
  {
    Cell_handle new_cell = tr.tds().create_cell();

    for (int v = 0; v < 4; v++){
      new_cell->set_vertex(v, fi.first->vertex(v));
    }

    new_cell->set_vertex(fi.second, vh);

    treat_new_cell(new_cell, subdomain, cell_selector, selected, c3t3);

    visitor.after_flip(new_cell);
    cells_to_update.push_back(new_cell);
  }

  //Update_existing cells
  for (const Facet& fi : facets_for_updated_cells)
  {
    fi.first->set_vertex(fi.second, vh);
    cells_to_update.push_back(fi.first);
  }

  typedef CGAL::Triple<Vertex_handle, Vertex_handle, Vertex_handle> Facet_vvv;
  typedef std::unordered_map<Facet_vvv, std::size_t> FaceMapIndex;

  FaceMapIndex facet_map_indices;
  std::vector<Facet> facets;

  for (const Facet& f : neighbor_facets)
  {
    Cell_handle ch = f.first;
    int v = f.second;

    Facet_vvv face = make_vertex_triple(ch->vertex(indices(v,0)),
                                        ch->vertex(indices(v,1)),
                                        ch->vertex(indices(v,2)));
    typename FaceMapIndex::iterator it = facet_map_indices.find(face);
    if (it == facet_map_indices.end())
    {
      facet_map_indices[face] = facets.size();
      facets.push_back(Facet(ch, v));
    }
  }

  //Update adjacencies and vertices cells
  for (Cell_handle ch : cells_to_update)
  {
    for (int v = 0; v < 4; v++)
    {
      Facet_vvv face = make_vertex_triple(ch->vertex(indices(v,0)),
                                          ch->vertex(indices(v,1)),
                                          ch->vertex(indices(v,2)));
      typename FaceMapIndex::iterator it = facet_map_indices.find(face);
      if (it == facet_map_indices.end())
      {
        facet_map_indices[face] = facets.size();
        facets.push_back(Facet(ch, v));
      }
      else
      {
        Facet facet = facets[it->second];

        //Update neighbor
        facet.first->set_neighbor(facet.second, ch);
        ch->set_neighbor(v, facet.first);
      }
      ch->vertex(v)->set_cell(ch);

      inc_cells[ch->vertex(v)].clear();
    }
    ch->reset_cache_validity();
  }

  // Update c3t3
  update_c3t3_facets(c3t3, cells_to_update, neighbor_facets);

  //Remove cells
  for (Cell_handle ch : to_remove)
  {
    treat_before_delete(ch, cell_selector, c3t3);
    ch->reset_cache_validity();
    tr.tds().delete_cell(ch);
  }

  ///********************VALIDITY CHECK***************************/
  //if (check_validity){

  //  double new_min_dh = DBL_MAX;

  //  post_sliver_Removal_cells.clear();
  //  for (unsigned int i = 0; i < cells_to_update.size(); ++i){
  //    post_sliver_Removal_cells.push_back(cells_to_update[i]);

  //    if (!tr.is_infinite(cells_to_update[i]))
  //      new_min_dh = std::min(new_min_dh, min_dihedral_angle(cells_to_update[i]));
  //  }

  //  post_sliver_Removal_vertices.clear();
  //  for (unsigned int i = 0; i < vertices_around_edge.size(); ++i){
  //    post_sliver_Removal_vertices.push_back(vertices_around_edge[i]);
  //  }

  //  current_edges.clear();
  //  for (unsigned int i = 0; i < vertices_around_edge.size(); ++i){
  //    current_edges.push_back(std::make_pair(vertices_around_edge[i]->point(), vh->point()));
  //  }


  //  for (unsigned int i = 0; i < cells_to_update.size(); ++i){
  //    if (!tr.is_valid(cells_to_update[i]))
  //      return INVALID_CELL;

  //    for (int v = 0; v < 4; v++){
  //      if (!tr.is_valid(cells_to_update[i]->neighbor(v)))
  //        return INVALID_CELL;

  //      if (!tr.tds().is_valid(cells_to_update[i]->vertex(v)))
  //        return INVALID_VERTEX;

  //    }
  //  }

  //  if ((current_min_dh - new_min_dh) > 0.01){
  //    std::cout << pre_sliver_Removal_cells.size() << " to " << post_sliver_Removal_cells.size() << " flip not improving the quality: " <<
  //      current_min_dh << " to " << new_min_dh << std::endl;
  //    return INVALID_CELL;
  //  }

  //}
  ///***********************************************************/

  // std::cout << "n_to_m_flip::end with success" << std::endl;

  return VALID_FLIP;
}


template<typename C3t3, typename IncCellsVectorMap, typename CellSelector, typename Visitor>
Sliver_removal_result flip_n_to_m(typename C3t3::Edge& edge,
                                  C3t3& c3t3,
                                  const std::vector<typename C3t3::Vertex_handle>& boundary_vertices,
                                  const Flip_Criterion& criterion,
                                  IncCellsVectorMap& inc_cells,
                                  CellSelector& cell_selector,
                                  Visitor& visitor)
{
  typedef typename C3t3::Vertex_handle Vertex_handle;
  typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
  typename C3t3::Triangulation& tr = c3t3.triangulation();

  Sliver_removal_result result = NOT_FLIPPABLE;

  typedef std::pair<Dihedral_angle_cosine, std::pair<Vertex_handle, int> > CosAngle_and_vertex;

  //std::cout << "n_to_m_flip " << boundary_vertices.size() << std::endl;
  if (criterion == MIN_ANGLE_BASED)
  {
    std::priority_queue<CosAngle_and_vertex,
                        std::vector<CosAngle_and_vertex>,
                        std::greater<CosAngle_and_vertex>
                      > candidates;

    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
    Cell_circulator done = circ;

    Dihedral_angle_cosine curr_max_cosdh = max_cos_dihedral_angle(tr, circ++);
    do
    {
      curr_max_cosdh = (std::max)(curr_max_cosdh, max_cos_dihedral_angle(tr, circ));
    } while (++circ != done);

    if (boundary_vertices.size() == 2)
      find_best_flip_to_improve_dh(c3t3, edge, boundary_vertices[0], boundary_vertices[1],
                                   candidates, curr_max_cosdh);
    else
      find_best_flip_to_improve_dh(c3t3, edge, candidates, curr_max_cosdh, inc_cells);

    bool flip_performed = false;
    while (!candidates.empty() && !flip_performed)
    {
      CosAngle_and_vertex curr_cost_vpair = candidates.top();
      candidates.pop();

//      std::cout << "\tcurrent   cos = " << curr_max_cosdh.value()
//        << "\t angle = " << std::acos(curr_max_cosdh.value()) * 180./CGAL_PI << std::endl;
//      std::cout << "\tcandidate cos = " << curr_cost_vpair.first.value()
//        << "\t angle = " << std::acos(curr_cost_vpair.first.value()) * 180./CGAL_PI << std::endl;
//      std::cout << std::endl;

      if (curr_max_cosdh <= curr_cost_vpair.first)
        return NO_BEST_CONFIGURATION;

      result = flip_n_to_m(c3t3, edge, curr_cost_vpair.second.first, inc_cells,
                           cell_selector, visitor);

      if (result != NOT_FLIPPABLE)
        flip_performed = true;
    }
  }

  return result;
}

template<typename C3t3, typename IncCellsVectorMap, typename CellSelector, typename Visitor>
Sliver_removal_result find_best_flip(typename C3t3::Edge& edge,
                                     C3t3& c3t3,
                                     const Flip_Criterion& criterion,
                                     IncCellsVectorMap& inc_cells,
                                     CellSelector& cell_selector,
                                     Visitor& visitor)
{
  typedef typename C3t3::Triangulation        Tr;
  typedef typename C3t3::Vertex_handle        Vertex_handle;
  typedef typename Tr::Facet_circulator       Facet_circulator;

  Tr& tr = c3t3.triangulation();

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

  Facet_circulator circ = tr.incident_facets(edge);
  Facet_circulator done = circ;

  //Identify the vertices around this edge
  std::unordered_set<Vertex_handle> vertices_around_edge;
  bool boundary_edge = false;
  bool hull_edge = false;

  std::unordered_set<Vertex_handle> boundary_vertices;
//  std::unordered_set<Vertex_handle> hull_vertices;
  do
  {
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle vi = circ->first->vertex(indices(circ->second, i));
      if (vi != v0 && vi != v1)
      {
        vertices_around_edge.insert(vi);

        if ( circ->first->subdomain_index()
             != circ->first->neighbor(circ->second)->subdomain_index())
        {
          boundary_edge = true;
          boundary_vertices.insert(vi);
        }

        if ( tr.is_infinite(circ->first)
             != tr.is_infinite(circ->first->neighbor(circ->second)))
        {
          hull_edge = true;
          //hull_vertices.insert(vi);
        }
      }
    }
  }
  while (++circ != done);


  //Check if not feature edge
  if (boundary_vertices.size() > 2)
    return NOT_FLIPPABLE;

  // perform flip when possible
  Sliver_removal_result res = NOT_FLIPPABLE;
  if (vertices_around_edge.size() == 3)
  {
    if (!boundary_edge && !hull_edge)
    {
      std::vector<Vertex_handle> vertices;
      vertices.insert(vertices.end(), vertices_around_edge.begin(), vertices_around_edge.end());
      res = flip_3_to_2(edge, c3t3, vertices, criterion, inc_cells, cell_selector);
    }
  }
  else
  {
    //TODO fix for hull edges
    // if( hull_edge )
    //    return n_to_m_flip( edge, hull_vertices, flip_criterion, check_validity );
    if (!hull_edge)
    {
      std::vector<Vertex_handle> vertices;
      vertices.insert(vertices.end(), boundary_vertices.begin(), boundary_vertices.end());
      res = flip_n_to_m(edge, c3t3, vertices, criterion, inc_cells, cell_selector, visitor);
      //return n_to_m_flip(edge, boundary_vertices, flip_criterion);
    }
  }

  return res;
}


template<typename VertexPair, typename C3t3,
         typename IncidentCellsVectorMap, typename CellSelector, typename Visitor>
std::size_t flip_all_edges(const std::vector<VertexPair>& edges,
                           C3t3& c3t3,
                           IncidentCellsVectorMap& inc_cells,
                           const Flip_Criterion& criterion,
                           CellSelector& cell_selector,
                           Visitor& visitor)
{
  typedef typename C3t3::Triangulation Tr;
  typedef typename Tr::Cell_handle   Cell_handle;
  typedef typename Tr::Edge          Edge;

  Tr& tr = c3t3.triangulation();

  std::size_t count = 0;
  for (const VertexPair& vp : edges)
  {
    boost::container::small_vector<Cell_handle, 64>& o_inc_vh = inc_cells[vp.first];
    if (o_inc_vh.empty())
    {
      CGAL_TR_LOCK_PROBE(tr, vp.first, "flip_edges_loop:vp.first");
      tr.incident_cells(vp.first, std::back_inserter(o_inc_vh));
    }

    Cell_handle ch;
    int i0, i1;
    if (is_edge_uv(vp.first, vp.second, o_inc_vh, ch, i0, i1))
    {
      Edge edge(ch, i0, i1);

      Sliver_removal_result res
        = find_best_flip(edge, c3t3, criterion, inc_cells, cell_selector, visitor);
      if (res == INVALID_CELL || res == INVALID_VERTEX || res == INVALID_ORIENTATION)
      {
        std::cout << "FLIP PROBLEM!!!!" << std::endl;
        return count;
      }
      if (res == VALID_FLIP)
      {
        ++count;
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
        std::cout << "\rFlip... (";
        std::cout << count << " flips)";
        std::cout.flush();
#endif
      }
    }
  }

  return count;
}

template<typename C3t3>
void collect_subdomains_on_boundary(const C3t3& c3t3,
  boost::unordered_map<typename C3t3::Vertex_handle,
    std::unordered_set<typename C3t3::Subdomain_index> >& vertices_subdomain_indices)
{
  for (auto c : c3t3.triangulation().all_cell_handles())
  {
    for (auto v : c3t3.triangulation().vertices(c))
    {
      const int dim = v->in_dimension();
      if(dim >= 0 && dim < 3)
        vertices_subdomain_indices[v].insert(c->subdomain_index());
    }
  }
}

template<typename C3T3, typename CellSelector>
void collectBoundaryEdgesAndComputeVerticesValences(
  const C3T3& c3t3,
  const CellSelector& cell_selector,
  std::vector<typename C3T3::Edge>& boundary_edges,
  boost::unordered_map<typename C3T3::Vertex_handle,
                       boost::unordered_map<typename C3T3::Surface_patch_index, unsigned int> >&
      boundary_vertices_valences,
  boost::unordered_map<typename C3T3::Vertex_handle, std::unordered_set<typename C3T3::Subdomain_index> >&
      vertices_subdomain_indices)
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef typename C3T3::Vertex_handle       Vertex_handle;
  typedef typename C3T3::Edge                Edge;
  typedef typename C3T3::Triangulation::Facet_circulator Facet_circulator;

  const typename C3T3::Triangulation& tr = c3t3.triangulation();

  boundary_edges.clear();
  boundary_vertices_valences.clear();

#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
  if (par_flip_collect())
  {
    boundary_edges = parallel_collect_finite_edges<Edge>(tr,
      [&c3t3, &cell_selector](const Edge& e, std::vector<Edge>& local)
      {
        if (is_boundary(c3t3, e, cell_selector))
          local.push_back(e);
      });
  }
  else
#endif
  for (const Edge& e : tr.finite_edges())
  {
    if (is_boundary(c3t3, e, cell_selector))
      boundary_edges.push_back(e);
  }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
  CGAL::Tetrahedral_remeshing::debug::dump_edges(boundary_edges,
                                                 "boundary_edges.polylines.txt");
#endif

  // collect incident subdomain indices at vertices
  collect_subdomains_on_boundary(c3t3, vertices_subdomain_indices);

  for (const Edge& e : boundary_edges)
  {
    const Vertex_handle v0 = e.first->vertex(e.second);
    const Vertex_handle v1 = e.first->vertex(e.third);

    const std::size_t n0 = vertices_subdomain_indices[v0].size();
    const std::size_t n1 = vertices_subdomain_indices[v1].size();

    //In case of feature edge
    if (n0 > 2 && n1 > 2)
    {
      Facet_circulator facet_circulator = tr.incident_facets(e);
      Facet_circulator done(facet_circulator);
      do
      {
        if (c3t3.is_in_complex(*facet_circulator))
        {
          Surface_patch_index surfi = c3t3.surface_patch_index(*facet_circulator);
          boundary_vertices_valences[v0][surfi]++;
          boundary_vertices_valences[v1][surfi]++;
        }
      } while (++facet_circulator != done);
    }
    //Normal surface edge, or non-manifold edge on dangling facet
    else
    {
      Facet_circulator facet_circulator = tr.incident_facets(e);
      Facet_circulator done(facet_circulator);
      Surface_patch_index first_patch = Surface_patch_index();
      do
      {
        if (c3t3.is_in_complex(*facet_circulator))
        {
          Surface_patch_index surfi = c3t3.surface_patch_index(*facet_circulator);
          if (first_patch == Surface_patch_index())
            first_patch = surfi;
          else if (first_patch == surfi)
            continue;

          boundary_vertices_valences[v0][surfi]++;
          boundary_vertices_valences[v1][surfi]++;
        }
      } while (++facet_circulator != done);
    }
  }
}

template<typename C3T3, typename IncCellsVector, typename Visitor>
Sliver_removal_result flip_n_to_m_on_surface(typename C3T3::Edge& edge,
    C3T3& c3t3,
    typename C3T3::Vertex_handle v0i,//v0 of new edge that will replace edge
    typename C3T3::Vertex_handle v1i,//v1 of new edge that will replace edge
    const IncCellsVector& cells_around_edge,
    Flip_Criterion /*flip_criterion*/,
    Visitor& /*visitor*/)
{
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef typename C3T3::Cell_handle   Cell_handle;

  typename C3T3::Triangulation& tr = c3t3.triangulation();

  const Vertex_handle u = edge.first->vertex(edge.second);
  const Vertex_handle v = edge.first->vertex(edge.third);

  typedef std::pair<int, int> IndInCell;
  std::map<Cell_handle, IndInCell> indices;
  for (Cell_handle c : cells_around_edge)
  {
    indices[c] = std::make_pair(c->index(u), c->index(v));
  }

  for (Cell_handle c : cells_around_edge)
  {
    int i = indices[c].first;
    int j = indices[c].second;
    c->set_vertex(i, v0i);
    c->set_vertex(j, v1i);

    if (!is_well_oriented(tr, c))
    {
      c->set_vertex(j, v0i);
      c->set_vertex(i, v1i);
      if (!is_well_oriented(tr, c))
      {
        //rollback all changes
        for (Cell_handle cc : cells_around_edge)
        {
          const int ii = indices[cc].first;
          if (cc->vertex(ii) != u)
          {
            cc->set_vertex(ii, u);
            cc->set_vertex(indices[cc].second, v);
          }
        }
        return NOT_FLIPPABLE;
      }
    }
  }

  for (Cell_handle c : cells_around_edge)
    c->reset_cache_validity();

  return VALID_FLIP;
}

//v0i and v1i are the vertices opposite to `edge`
//on facets of the surface
template<typename C3T3, typename IncCellsVectorMap,
         typename Visitor>
Sliver_removal_result flip_on_surface(C3T3& c3t3,
    typename C3T3::Edge& edge,
    const typename C3T3::Vertex_handle v0i,//v0 of new edge that will replace edge
    const typename C3T3::Vertex_handle v1i,//v1 of new edge that will replace edge
    IncCellsVectorMap& inc_cells,
    Flip_Criterion flip_criterion,
    Visitor& visitor)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename Tr::Cell_handle     Cell_handle;
  typedef typename Tr::Vertex_handle   Vertex_handle;
  typedef typename Tr::Cell_circulator Cell_circulator;

  Tr& tr = c3t3.triangulation();
  Cell_circulator circ = tr.incident_cells(edge);
  Cell_circulator done(circ);

  std::vector<Cell_handle> cells_around_edge;
  do
  {
    cells_around_edge.push_back(circ);
  } while (++circ != done);

  if (cells_around_edge.size() != 4)
  {
    if (cells_around_edge.size() > 4){
//////      if (flip_criterion == VALENCE_BASED){
//////        return find_best_n_m_flip(edge, vh0_index, vh1_index);
//////      }
//////      else {
//        std::vector<Vertex_handle> boundary_vertices;
//        boundary_vertices.push_back(v0i);
//        boundary_vertices.push_back(v1i);
#ifdef CGAL_FLIP_ON_SURFACE_DISABLE_NM_FLIP
        return NOT_FLIPPABLE;
#else
        return flip_n_to_m_on_surface(edge, c3t3, v0i, v1i,
                                      cells_around_edge, flip_criterion,
                                      visitor);
#endif
//////    }
    }
    else
      return NOT_FLIPPABLE;
  }

#ifdef CGAL_FLIP_ON_SURFACE_DISABLE_44_FLIP
  return NOT_FLIPPABLE;
#endif

  inc_cells[edge.first->vertex(edge.second)].clear();
  inc_cells[edge.first->vertex(edge.third)].clear();

  Cell_handle ch0, ch1, ch2, ch3;
  ch0 = cells_around_edge[0];
  ch1 = cells_around_edge[1];
  ch2 = cells_around_edge[2];
  ch3 = cells_around_edge[3];

  Dihedral_angle_cosine curr_max_cosdh = max_cos_dihedral_angle(tr, ch0);
  for (int i = 1; i < 4; ++i)
    curr_max_cosdh = (std::max)(curr_max_cosdh,
                                max_cos_dihedral_angle(tr, cells_around_edge[i]));

  Vertex_handle vh0, vh1, vh2, vh3, vh4, vh5;

  int ivh4_in_ch0 = ch0->index(ch1);
  vh4 = ch0->vertex(ivh4_in_ch0);

  int ivh2_in_ch0 = ch0->index(ch3);
  vh2 = ch0->vertex(ivh2_in_ch0);

  vh5 = ch1->vertex(ch1->index(ch0));
  vh0 = ch2->vertex(ch2->index(ch1));

  for (int j = 0; j < 3; j++){
    if (indices(ivh4_in_ch0, j) == ivh2_in_ch0){
      int j1 = (j + 1) % 3;
      int j2 = (j + 2) % 3;
      vh1 = ch0->vertex(indices(ivh4_in_ch0, j1));
      vh3 = ch0->vertex(indices(ivh4_in_ch0, j2));
      break;
    }
  }

  bool planar_flip;
  if ((vh0 == v0i && vh2 == v1i) || (vh2 == v0i && vh0 == v1i))
    planar_flip = true;
  else if ((vh4 == v0i && vh5 == v1i) || (vh5 == v0i && vh4 == v1i))
    planar_flip = false;
  else
    return NOT_FLIPPABLE;

  typedef typename C3T3::Facet Facet;
  typedef typename C3T3::Surface_patch_index Surface_patch_index;

  if (planar_flip)
  {
#ifdef CGAL_FLIP_ON_SURFACE_DISABLE_PLANAR_44_FLIP
    return NOT_FLIPPABLE;
#endif
    Surface_patch_index patch = c3t3.surface_patch_index(ch0, ch0->index(vh4));
    CGAL_assertion(patch != Surface_patch_index());
    CGAL_assertion(c3t3.is_in_complex(ch0, ch0->index(vh4)));
    c3t3.remove_from_complex(ch0, ch0->index(vh4));
    CGAL_assertion(c3t3.is_in_complex(ch3, ch3->index(vh4)));
    c3t3.remove_from_complex(ch3, ch3->index(vh4));

    boost::unordered_map<Facet, Surface_patch_index> opposite_facet_in_complex;
    for (Cell_handle chi : cells_around_edge)
    {
      Facet f1(chi, chi->index(vh1));
      Facet f2(chi, chi->index(vh3));

      if (c3t3.is_in_complex(f1))
      {
        Surface_patch_index spi = c3t3.surface_patch_index(f1);
        opposite_facet_in_complex[c3t3.triangulation().mirror_facet(f1)] = spi;
        c3t3.remove_from_complex(f1);
      }
      if (c3t3.is_in_complex(f2))
      {
        Surface_patch_index spi = c3t3.surface_patch_index(f2);
        opposite_facet_in_complex[c3t3.triangulation().mirror_facet(f2)] = spi;
        c3t3.remove_from_complex(f2);
      }
    }

    Cell_handle n_ch3_vh1 = ch3->neighbor(ch3->index(vh1));
    Cell_handle n_ch0_vh3 = ch0->neighbor(ch0->index(vh3));

    Cell_handle n_ch2_vh1 = ch2->neighbor(ch2->index(vh1));
    Cell_handle n_ch1_vh3 = ch1->neighbor(ch1->index(vh3));

    ch3->set_vertex(ch3->index(vh3), vh2);
    ch0->set_vertex(ch0->index(vh1), vh0);
    ch2->set_vertex(ch2->index(vh3), vh2);
    ch1->set_vertex(ch1->index(vh1), vh0);

    Sliver_removal_result db = VALID_FLIP;
    if (!is_well_oriented(tr, ch0)
      || !is_well_oriented(tr, ch1)
      || !is_well_oriented(tr, ch2)
      || !is_well_oriented(tr, ch3))
      db = NOT_FLIPPABLE;
    else if (curr_max_cosdh < max_cos_dihedral_angle(tr, ch0, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch1, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch2, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch3, false))
      db = NO_BEST_CONFIGURATION;

    if(db != VALID_FLIP)
    {
      ch3->set_vertex(ch3->index(vh2), vh3);
      ch0->set_vertex(ch0->index(vh0), vh1);
      ch2->set_vertex(ch2->index(vh2), vh3);
      ch1->set_vertex(ch1->index(vh0), vh1);

      c3t3.add_to_complex(ch0, ch0->index(vh4), patch);
      c3t3.add_to_complex(ch3, ch3->index(vh4), patch);

      for (Cell_handle chi : cells_around_edge)
      {
        Facet f1(chi, chi->index(vh1));
        Facet f2(chi, chi->index(vh3));

        auto it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f1));
        if (it != opposite_facet_in_complex.end())
          c3t3.add_to_complex(f1, it->second);

        it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f2));
        if (it != opposite_facet_in_complex.end())
          c3t3.add_to_complex(f2, it->second);
      }

      return db;
    }

    //Top cells 2-2 flip
    ch3->set_neighbor(ch3->index(vh1), ch0);
    ch3->set_neighbor(ch3->index(vh0), n_ch0_vh3);
    n_ch0_vh3->set_neighbor(n_ch0_vh3->index(ch0), ch3);

    ch0->set_neighbor(ch0->index(vh3), ch3);
    ch0->set_neighbor(ch0->index(vh2), n_ch3_vh1);
    n_ch3_vh1->set_neighbor(n_ch3_vh1->index(ch3), ch0);

    //Bottom cells 2-2 flip
    ch2->set_neighbor(ch2->index(vh1), ch1);
    ch2->set_neighbor(ch2->index(vh0), n_ch1_vh3);
    n_ch1_vh3->set_neighbor(n_ch1_vh3->index(ch1), ch2);

    ch1->set_neighbor(ch1->index(vh3), ch2);
    ch1->set_neighbor(ch1->index(vh2), n_ch2_vh1);
    n_ch2_vh1->set_neighbor(n_ch2_vh1->index(ch2), ch1);

    for (Cell_handle ci : cells_around_edge)
    {
      for (int j = 0; j < 4; j++)
        ci->vertex(j)->set_cell(ci);
    }

    c3t3.add_to_complex(ch0, ch0->index(vh4), patch);
    c3t3.add_to_complex(ch3, ch3->index(vh4), patch);

    for (Cell_handle chi : cells_around_edge)
    {
      Facet f1(chi, chi->index(vh0));
      Facet f2(chi, chi->index(vh2));

      auto it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f1));
      if (it != opposite_facet_in_complex.end())
        c3t3.add_to_complex(f1, it->second);

      it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f2));
      if (it != opposite_facet_in_complex.end())
        c3t3.add_to_complex(f2, it->second);
    }

    for(Cell_handle c : cells_around_edge)
      c->reset_cache_validity();

    return db;
  }
  else //Non planar flip
  {
#ifdef CGAL_FLIP_ON_SURFACE_DISABLE_NON_PLANAR_44_FLIP
    return NOT_FLIPPABLE;
#endif
    typename C3T3::Surface_patch_index patch = c3t3.surface_patch_index(ch0, ch0->index(vh2));
    CGAL_assertion(patch != typename C3T3::Surface_patch_index());

    CGAL_assertion(c3t3.is_in_complex(ch0, ch0->index(vh2)));
    c3t3.remove_from_complex(ch0, ch0->index(vh2));
    CGAL_assertion(c3t3.is_in_complex(ch1, ch1->index(vh2)));
    c3t3.remove_from_complex(ch1, ch1->index(vh2));

    boost::unordered_map<Facet, Surface_patch_index> opposite_facet_in_complex;
    for (Cell_handle chi : cells_around_edge)
    {
      Facet f1(chi, chi->index(vh1));
      Facet f2(chi, chi->index(vh3));

      if (c3t3.is_in_complex(f1))
      {
        Surface_patch_index spi = c3t3.surface_patch_index(f1);
        opposite_facet_in_complex[c3t3.triangulation().mirror_facet(f1)] = spi;
        c3t3.remove_from_complex(f1);
      }
      if (c3t3.is_in_complex(f2))
      {
        Surface_patch_index spi = c3t3.surface_patch_index(f2);
        opposite_facet_in_complex[c3t3.triangulation().mirror_facet(f2)] = spi;
        c3t3.remove_from_complex(f2);
      }
    }

    // Top Flip
    ch3->set_vertex(ch3->index(vh1), vh5);
    ch2->set_vertex(ch2->index(vh3), vh4);
    ch0->set_vertex(ch0->index(vh1), vh5);
    ch1->set_vertex(ch1->index(vh3), vh4);

    Sliver_removal_result db = VALID_FLIP;
    if (!is_well_oriented(tr, ch0)
      || !is_well_oriented(tr, ch1)
      || !is_well_oriented(tr, ch2)
      || !is_well_oriented(tr, ch3))
      db = NOT_FLIPPABLE;
    else if (curr_max_cosdh < max_cos_dihedral_angle(tr, ch0, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch1, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch2, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch3, false))
      db = NO_BEST_CONFIGURATION;

    if (db == NOT_FLIPPABLE || db == NO_BEST_CONFIGURATION)
    {
      ch3->set_vertex(ch3->index(vh5), vh1);
      ch2->set_vertex(ch2->index(vh4), vh3);
      ch0->set_vertex(ch0->index(vh5), vh1);
      ch1->set_vertex(ch1->index(vh4), vh3);

      c3t3.add_to_complex(ch0, ch0->index(vh2), patch);
      c3t3.add_to_complex(ch1, ch1->index(vh2), patch);

      for (Cell_handle chi : cells_around_edge)
      {
        Facet f1(chi, chi->index(vh1));
        Facet f2(chi, chi->index(vh3));

        auto it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f1));
        if (it != opposite_facet_in_complex.end())
          c3t3.add_to_complex(f1, it->second);

        it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f2));
        if (it != opposite_facet_in_complex.end())
          c3t3.add_to_complex(f2, it->second);
      }

      return db;
    }

    //Left cells 2-2 flip
    Cell_handle n_ch3_vh3 = ch3->neighbor(ch3->index(vh3));
    Cell_handle n_ch2_vh1 = ch2->neighbor(ch2->index(vh1));

    ch3->set_neighbor(ch3->index(vh3), ch2);
    ch3->set_neighbor(ch3->index(vh4), n_ch2_vh1);
    n_ch2_vh1->set_neighbor(n_ch2_vh1->index(ch2), ch3);

    ch2->set_neighbor(ch2->index(vh1), ch3);
    ch2->set_neighbor(ch2->index(vh5), n_ch3_vh3);
    n_ch3_vh3->set_neighbor(n_ch3_vh3->index(ch3), ch2);

    //Right cells 2-2 flip
    Cell_handle n_ch0_vh3 = ch0->neighbor(ch0->index(vh3));
    Cell_handle n_ch1_vh1 = ch1->neighbor(ch1->index(vh1));

    ch0->set_neighbor(ch0->index(vh3), ch1);
    ch0->set_neighbor(ch0->index(vh4), n_ch1_vh1);
    n_ch1_vh1->set_neighbor(n_ch1_vh1->index(ch1), ch0);

    ch1->set_neighbor(ch1->index(vh1), ch0);
    ch1->set_neighbor(ch1->index(vh5), n_ch0_vh3);
    n_ch0_vh3->set_neighbor(n_ch0_vh3->index(ch0), ch1);

    for (const Cell_handle ce : cells_around_edge)
    {
      for (int j = 0; j < 4; j++)
        ce->vertex(j)->set_cell(ce);
    }

    c3t3.add_to_complex(ch0, ch0->index(vh2), patch);
    c3t3.add_to_complex(ch1, ch1->index(vh2), patch);

    for (Cell_handle chi : cells_around_edge)
    {
      Facet f1(chi, chi->index(vh4));
      Facet f2(chi, chi->index(vh5));

      auto it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f1));
      if (it != opposite_facet_in_complex.end())
        c3t3.add_to_complex(f1, it->second);

      it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f2));
      if (it != opposite_facet_in_complex.end())
        c3t3.add_to_complex(f2, it->second);
    }

    for (Cell_handle c : cells_around_edge)
      c->reset_cache_validity();

    return VALID_FLIP;
  }

  return NOT_FLIPPABLE;
}

template<typename C3T3, typename SurfaceIndexMapMap,
         typename IncidentCellsVectorMap, typename Flip_Criterion,
         typename CellSelector, typename Visitor>
std::size_t flipBoundaryEdges(C3T3& c3t3,
                              const std::vector<typename C3T3::Edge>& boundary_edges,
                              SurfaceIndexMapMap& boundary_vertices_valences,
                              IncidentCellsVectorMap& inc_cells,
                              const Flip_Criterion& flip_criterion,
                              CellSelector& cell_selector,
                              Visitor& visitor)
{
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef typename C3T3::Facet Facet;
  typedef typename C3T3::Edge Edge;
  typedef typename C3T3::Triangulation Tr;
  typedef std::pair<Vertex_handle, Vertex_handle> Edge_vv;

  std::size_t nb_success = 0;

  Tr& tr = c3t3.triangulation();

  std::vector<Edge_vv> candidate_edges_for_flip;
  for(const Edge& e : boundary_edges) {
    if(!c3t3.is_in_complex(e))
      candidate_edges_for_flip.push_back(make_vertex_pair(e));
  }

  for(const auto& [vh0, vh1] : candidate_edges_for_flip) {
    boost::container::small_vector<Cell_handle, 64>& inc_vh0 = inc_cells[vh0];
    if(inc_vh0.empty())
    {
      CGAL_TR_LOCK_PROBE(tr, vh0, "flip_all_edges:vh0");
      tr.incident_cells(vh0, std::back_inserter(inc_vh0));
    }

    Cell_handle c;
    int i, j;
    if(!is_edge_uv(vh0, vh1, inc_vh0, c, i, j))
      continue;

    Edge edge(c, i, j);
    std::vector<Facet> boundary_facets;
    const bool on_boundary = is_boundary_edge(edge, c3t3, cell_selector, boundary_facets);

    //    if (on_boundary && boundary_facets.empty())
    //    {
    //      std::cerr << vh0->point().point() << "\t " << vh1->point().point() << std::endl;
    //      bool b = is_boundary_edge(vh0, vh1, c3t3, cell_selector);
    //      CGAL::Tetrahedral_remeshing::debug::dump_c3t3(c3t3, "dump_c3t3_about_boundary_");
    //      CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(c3t3, "dump_facets_about_boundary_.off");
    //      CGAL::Tetrahedral_remeshing::debug::dump_facets_from_selection(
    //        c3t3, cell_selector, "dump_facets_from_selection_.off");
    //      std::cerr << "valid = " << tr.tds().is_valid(true) << std::endl;
    //      std::cerr << "boundary = " << b << std::endl;
    //      CGAL_assertion(on_boundary);
    //    }
    //    else if (on_boundary && boundary_facets.size() != 2)
    //    {
    //      std::cerr << vh0->point().point() << "\t " << vh1->point().point() << std::endl;
    //      CGAL::Tetrahedral_remeshing::debug::dump_c3t3(c3t3, "dump_c3t3_about_boundary_");
    //      CGAL::Tetrahedral_remeshing::debug::dump_facets(boundary_facets, "dump_boundary_facets.polylines.txt");
    //      std::vector<Facet> dummy_facets;
    //      bool b = is_boundary_edge(edge, c3t3, cell_selector, dummy_facets, true/**/);
    //      std::cerr << "boundary = " << b << std::endl;
    //    }

    if(!on_boundary || boundary_facets.size() != 2)
      continue;
    CGAL_assertion(boundary_facets.size() == 2);

    if(flip_surface_edge(c3t3, edge, boundary_facets, boundary_vertices_valences,
                          inc_cells, flip_criterion, visitor))
    {
      ++nb_success;
      CGAL_expensive_assertion(tr.tds().is_valid());
    }
  }
  return nb_success;
}

//todo : write this function
template <typename C3t3,
          typename BV_valences,
          typename IncidentCellsVectorMap,
          typename Flip_criterion,
          typename Visitor>
bool flip_surface_edge(C3t3& c3t3,
                        typename C3t3::Edge& edge,
                        const std::vector<typename C3t3::Facet>& boundary_facets,
                        BV_valences& boundary_vertices_valences,
                        IncidentCellsVectorMap& inc_cells,
                        const Flip_criterion& flip_criterion,
                        Visitor& visitor)
{
  using Vertex_handle       = typename C3t3::Vertex_handle;
  using Facet               = typename C3t3::Facet;
  using Cell_handle         = typename C3t3::Cell_handle;
  using Surface_patch_index = typename C3t3::Surface_patch_index;

  auto& tr = c3t3.triangulation();

    const Facet& f0 = boundary_facets[0];
    const Facet& f1 = boundary_facets[1];

    const Vertex_handle vh0 = edge.first->vertex(edge.second);
    const Vertex_handle vh1 = edge.first->vertex(edge.third);

    // find 3rd and 4th vertices to flip on surface
    const Vertex_handle vh2 = third_vertex(f0, vh0, vh1, tr);
    const Vertex_handle vh3 = third_vertex(f1, vh0, vh1, tr);

    CGAL_expensive_assertion(debug::check_facets(vh0, vh1, vh2, vh3, c3t3));

    if (!tr.tds().is_edge(vh2, vh3)) // most-likely to happen early exit
    {
      const Surface_patch_index surfi = c3t3.surface_patch_index(boundary_facets[0]);

      int v0 = boundary_vertices_valences.at(vh0)[surfi];
      int v1 = boundary_vertices_valences.at(vh1)[surfi];
      int v2 = boundary_vertices_valences.at(vh2)[surfi];
      int v3 = boundary_vertices_valences.at(vh3)[surfi];

      if(v0 < 2 || v1 < 2 || v2 < 2 || v3 < 2)
        return false;

      int m0 = (boundary_vertices_valences.at(vh0).size() > 1 ? 4 : 6);
      int m1 = (boundary_vertices_valences.at(vh1).size() > 1 ? 4 : 6);
      int m2 = (boundary_vertices_valences.at(vh2).size() > 1 ? 4 : 6);
      int m3 = (boundary_vertices_valences.at(vh3).size() > 1 ? 4 : 6);

      int initial_cost = (v0 - m0)*(v0 - m0)
                       + (v1 - m1)*(v1 - m1)
                       + (v2 - m2)*(v2 - m2)
                       + (v3 - m3)*(v3 - m3);
      v0--;
      v1--;
      v2++;
      v3++;

      int final_cost = (v0 - m0)*(v0 - m0)
                     + (v1 - m1)*(v1 - m1)
                     + (v2 - m2)*(v2 - m2)
                     + (v3 - m3)*(v3 - m3);
      if (initial_cost > final_cost)
      {
        CGAL_expensive_assertion_code(std::size_t nbf =
          std::distance(c3t3.facets_in_complex_begin(),
                        c3t3.facets_in_complex_end()));
        CGAL_expensive_assertion_code(std::size_t nbe =
          std::distance(c3t3.edges_in_complex_begin(),
                        c3t3.edges_in_complex_end()));

        Sliver_removal_result db = flip_on_surface(c3t3, edge, vh2, vh3,
                                                   inc_cells,
                                                   flip_criterion,
                                                   visitor);
        if (db == VALID_FLIP)
        {
          CGAL_expensive_assertion(tr.tds().is_edge(vh2, vh3));
          Cell_handle c;
          int li, lj, lk;
          CGAL_expensive_assertion_code(bool b =)
          tr.tds().is_facet(vh2, vh3, vh0, c, li, lj, lk);
          CGAL_expensive_assertion(b);
          c3t3.add_to_complex(c, (6 - li - lj - lk), surfi);

          CGAL_expensive_assertion_code(b = )
          tr.tds().is_facet(vh2, vh3, vh1, c, li, lj, lk);
          CGAL_expensive_assertion(b);
          c3t3.add_to_complex(c, (6 - li - lj - lk), surfi);

          CGAL_expensive_assertion_code(std::size_t nbf_post =
            std::distance(c3t3.facets_in_complex_begin(),
                          c3t3.facets_in_complex_end()));
          CGAL_expensive_assertion(nbf == nbf_post);
          CGAL_expensive_assertion_code(std::size_t nbe_post =
            std::distance(c3t3.edges_in_complex_begin(),
                          c3t3.edges_in_complex_end()));
          CGAL_expensive_assertion(nbe == nbe_post);

          boundary_vertices_valences[vh0][surfi]--;
          boundary_vertices_valences[vh1][surfi]--;
          boundary_vertices_valences[vh2][surfi]++;
          boundary_vertices_valences[vh3][surfi]++;

          return true;
        }
      }
    }
  return false;
}

// Shared state for the internal and boundary edge-flip operations: the cell
// selector, the visitor, and the incident-cells cache used by find_best_flip
// and flip_on_surface.
template <typename C3t3,
          typename CellSelector,
          typename Visitor>
class Edge_flip_operation_base
{
protected:
  using Tr            = typename C3t3::Triangulation;
  using Cell_handle   = typename Tr::Cell_handle;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Edge          = typename Tr::Edge;
  using Facet         = typename Tr::Facet;
  using Cells_vector  = boost::container::small_vector<Cell_handle, 64>;
  // Per-worker incident-cell cache. It MUST be per-worker rather than one shared
  // concurrent map: lock-elided interior flips run execute_operation lock-free,
  // and find_best_flip lazily caches inc_cells[w] for the opposite/ring vertices w
  // around the edge -- vertices the interior/boundary classification neither tags
  // nor locks. Two elided flips in adjacent buckets that share such a w would
  // otherwise concurrently mutate the same Cells_vector (a small_vector, not
  // thread-safe) -> heap corruption -> rare, delayed segfault or a corrupted
  // triangulation that only surfaces phases later. lock_zone (populate) and
  // execute_operation (read/cache) run on the same worker, so a thread-local map
  // keeps the whole cache benefit with no cross-thread sharing.
  using Incident_cells_local = std::unordered_map<Vertex_handle, Cells_vector>;
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
  using Incident_cells_map = tbb::enumerable_thread_specific<Incident_cells_local>;
#else
  using Incident_cells_map = Incident_cells_local;
#endif

  CellSelector& m_cell_selector;
  Visitor& m_visitor;
  // Shared across the internal and boundary flip passes, exactly as the former
  // flip_edges() shared a single inc_cells map between flip_all_edges() and
  // flipBoundaryEdges().
  Incident_cells_map& inc_cells;

  Edge_flip_operation_base(CellSelector& cell_selector,
                        Visitor& visitor,
                        Incident_cells_map& incident_cells)
      : m_cell_selector(cell_selector)
      , m_visitor(visitor)
      , inc_cells(incident_cells) {}

  // This worker's slice of the shared cache (see Incident_cells_map above).
  Incident_cells_local& inc_cells_map() const
  {
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
    return inc_cells.local();
#else
    return inc_cells;
#endif
  }

public:
  typename Tr::Geom_traits::Point_3
  point_on_element(const std::pair<Vertex_handle, Vertex_handle>& vp) const
  {
    auto cp = typename Tr::Geom_traits().construct_point_3_object();
    return cp(vp.first->point());
  }
};

// Flip of internal (non-boundary) edges. Mirrors the former flip_all_edges():
// reset the cell caches, collect the internal edges, then run find_best_flip
// on each in turn, sharing the incident-cells cache.
template <typename C3t3, typename CellSelector, typename Visitor>
class Internal_edge_flip_operation
    : public Edge_flip_operation_base<C3t3, CellSelector, Visitor>,
      public Elementary_operation<C3t3,
                                 std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                                 std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>>
{
  using BaseClass = Edge_flip_operation_base<C3t3, CellSelector, Visitor>;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::Cells_vector;
  using BaseClass::m_cell_selector;
  using BaseClass::m_visitor;
  using BaseClass::inc_cells;
  // Dependent base member: GCC/Clang need the using-declaration, MSVC does not.
  using BaseClass::inc_cells_map;

public:
  using Incident_cells_map = typename BaseClass::Incident_cells_map;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using Base_operation = Elementary_operation<C3t3, Edge_vv, std::vector<Edge_vv>>;
  using Element_type = typename Base_operation::Element_type;
  static_assert(std::is_same_v<Element_type, Edge_vv>, "Element_type must be Edge_vv");
  using Element_range = typename Base_operation::Element_range;

  Internal_edge_flip_operation(CellSelector& cell_selector,
                            Visitor& visitor,
                            Incident_cells_map& incident_cells)
      : BaseClass(cell_selector, visitor, incident_cells) {}

  Element_range get_elements(const C3t3& c3t3) const override
  {
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
    if (par_flip_collect())
    {
      const auto& tr = c3t3.triangulation();

      // we will use sliver_value to store the cos_dihedral_angle
      std::vector<typename C3t3::Cell_handle> complex_cells;
      complex_cells.reserve(c3t3.number_of_cells_in_complex() + 64);
      for (auto c : c3t3.cells_in_complex())
        complex_cells.push_back(c);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, complex_cells.size()),
        [&complex_cells](const tbb::blocked_range<std::size_t>& r)
        {
          for (std::size_t i = r.begin(); i != r.end(); ++i)
            complex_cells[i]->reset_cache_validity();
        });

      const CellSelector& sel = m_cell_selector;
      return parallel_collect_finite_edges<Edge_vv>(tr,
        [&c3t3, &sel](const typename C3t3::Triangulation::Edge& e,
                      std::vector<Edge_vv>& local)
        {
          if (is_internal(e, c3t3, sel))
            local.push_back(make_vertex_pair(e));
        });
    }
#endif

    for (auto c : c3t3.cells_in_complex())
      c->reset_cache_validity();//we will use sliver_value
                                //to store the cos_dihedral_angle

    std::vector<Edge_vv> inside_edges;
    get_internal_edges(c3t3, m_cell_selector, std::back_inserter(inside_edges));
    return inside_edges;
  }

  bool execute_operation(const Element_type& vp, C3t3& c3t3) override
  {
    // Any read of a vertex star that cannot acquire its lock throws
    // Flip_lock_bail. All such reads happen during selection, before the
    // triangulation is modified, so abandoning here is always safe: the edge
    // is simply left unflipped for a later pass to retry.
    try
    {
      auto& tr = c3t3.triangulation();

      // Both endpoints are read and modified by every flip on this edge.
      CGAL_TR_REQUIRE_LOCK(tr, vp.first);
      CGAL_TR_REQUIRE_LOCK(tr, vp.second);

      // The per-worker cache only stays valid while this worker holds the
      // locks under which it was built. unlock_all_elements() after the
      // previous element dropped those, so any surviving entry may name a
      // cell another worker has since deleted. Start each element clean.
      // No flush needed per element: locks are held for the whole bucket
      // (prefers_bucket_scoped_locks()), so entries cached under them stay
      // valid until on_locks_released() drops them.
      auto& inc_cells_l = inc_cells_map();

      Cells_vector& o_inc_vh = inc_cells_l[vp.first];
      if (o_inc_vh.empty())
      {
        CGAL_TR_LOCK_PROBE(tr, vp.first,
                           "Internal_flip::execute_operation:vp.first");
        tr.incident_cells(vp.first, std::back_inserter(o_inc_vh));
      }

      Cell_handle ch;
      int i0, i1;
      if (!is_edge_uv(vp.first, vp.second, o_inc_vh, ch, i0, i1))
        return false;

      Edge edge(ch, i0, i1);
      const Sliver_removal_result res
        = find_best_flip(edge, c3t3, MIN_ANGLE_BASED, inc_cells_l, m_cell_selector, m_visitor);
      return (res == VALID_FLIP);
    }
    catch (const Flip_lock_bail&)
    {
      // Contended neighbourhood: give up on this edge, keep the mesh valid.
      // Release what we hold so the contending worker can proceed; that
      // invalidates everything cached under those locks.
      CGAL_TR_COUNT(lock_retry);
      c3t3.triangulation().unlock_all_elements();
      on_locks_released();
      return false;
    }
  }

  // Not used any more: prefers_conflict_free_waves() below routes this
  // operation through apply_conflict_free_waves() instead of
  // apply_unordered_processing(), so lock_zone() is never called. Kept in its
  // original, simple, safe form rather than deleted -- two more elaborate
  // versions were tried and rejected here earlier this session: locking only
  // {e.first, e.second} left find_best_flip_to_improve_dh()/flip_n_to_m()'s
  // ring-vertex reads unprotected regardless of classification (confirmed by
  // a SIGSEGV backtrace inside find_best_flip_to_improve_dh on
  // mesh3_243015); widening it to also lock ring vertices closed that but
  // required either an unordered or an ordered multi-vertex acquisition, and
  // both ran into the same problem: safely discovering *which* ring vertices
  // to lock requires reading the triangulation around the edge, which is
  // itself unsafe before any lock is held (confirmed by a SIGSEGV inside
  // Triangulation_data_structure_3::is_edge(), called from inside a from-
  // scratch lock_zone() attempt with no lock held yet). See the
  // "mesh3_243015 flip defect" log entries (third and fourth occurrence) for
  // both attempts and why apply_conflict_free_waves() sidesteps the
  // bootstrapping problem entirely: its claim is computed serially during
  // wave selection, before any parallel mutation starts, so there is nothing
  // to race against yet.
  bool lock_zone(const Element_type& e, const C3t3& c3t3) const override
  {
    auto& tr = c3t3.triangulation();
    Cells_vector inc_cells_first, inc_cells_second;
    const bool locked = tr.try_lock_and_get_incident_cells(e.first, inc_cells_first)
                      && tr.try_lock_and_get_incident_cells(e.second, inc_cells_second);
    if (locked)
    {
      auto& inc_cells_l = inc_cells_map(); // per-worker cache
      inc_cells_l[e.first] = inc_cells_first;
      inc_cells_l[e.second] = inc_cells_second;
    }
    return locked;
  }

  // *** OPEN DEFECT, NOT CLOSED after SIX attempts -- see the "mesh3_243015
  // flip defect" log entries for the full record. Restored (again) to the
  // pre-session form: {e.first, e.second} only, no lock_zone_precomputed
  // override (falls back to lock_zone() above, unchanged). Flip's true
  // footprint DOES reach its edge's ring (opposite) vertices --
  // find_best_flip_to_improve_dh() walks incident_cells(w) for every ring
  // vertex w, and flip_n_to_m() rewrites cells incident to w (this file) --
  // and this under-claims relative to that. What was tried and rejected,
  // in order:
  // 1. Widen this method to include ring vertices, keep
  //    apply_unordered_processing()'s lock_zone() at {e.first, e.second}
  //    only: classification correctly routes affected candidates to the
  //    locked path, but lock_zone() itself only ever locks the edge's own
  //    endpoints, so two "correctly classified" candidates sharing a ring
  //    vertex can still race on it. Confirmed live: a SIGSEGV backtrace
  //    inside find_best_flip_to_improve_dh().
  // 2. Widen lock_zone() to also lock ring vertices, unordered (v0, v1, then
  //    ring vertices in circulation order): genuine multi-GB memory growth
  //    under contention (no consistent lock order across threads -> two
  //    candidates can each hold what the other wants, both fail, unlock,
  //    retry -- with every failed attempt still paying for each vertex's
  //    incident-cell fetch before discovering the conflict).
  // 3. Same, but with a fixed global order (sort the vertex set by &*v
  //    before acquiring): the ordering fix is sound in principle, but
  //    *discovering* the ring-vertex set itself requires reading the
  //    triangulation around the edge, which is unsafe before any lock is
  //    held. Confirmed live: a SIGSEGV inside
  //    Triangulation_data_structure_3::is_edge(), called from lock_zone()
  //    with no lock held yet.
  // 4. Route through apply_conflict_free_waves() instead (its claim is
  //    computed serially during wave selection, sidestepping lock_zone()'s
  //    bootstrapping problem entirely, and is provably sound given a widened
  //    locked_vertices()): correct, but impractically slow. Internal flip's
  //    candidate set is essentially every internal edge (dense, unlike
  //    collapse's short-edge subset), so claims overlap constantly and very
  //    few candidates are accepted per wave -- over 2.5 CPU-minutes at
  //    ~100% (wave selection is serial) without finishing one iteration's
  //    flip phase on mesh3_243015, against a normal 10-40s run.
  // 5. Widen only this method (classification), leave
  //    apply_unordered_processing()'s original single-layer
  //    footprint_needs_halo() growth as-is (attempt 1's premise, but
  //    actually tested this time -- attempt 1 above was measured against a
  //    version of this method that had been accidentally reverted mid-session
  //    by an edit meant for Boundary_edge_flip_operation, so it never really
  //    ran): still crashed with std::bad_alloc, preceded by a memory spike
  //    (~10GB) far too large to be legitimate bookkeeping for a ~30K-vertex
  //    mesh -- consistent with the same class of corruption as before,
  //    meaning even the narrowest version of this widening does not close
  //    the gap on its own.
  //
  // 6. Route Pass 3's boundary path through a NEW lock_zone_precomputed()
  //    hook instead of lock_zone(), given the exact vertex small_vector Pass 1
  //    already computed (safely, mesh read-only) via a widened
  //    locked_vertices() -- so Pass 3 never re-discovers the ring set live,
  //    only locks handles it already knows are valid. This does sidestep
  //    attempt 3's SIGSEGV (no live is_edge() call at Pass-3 time). It still
  //    hit std::bad_alloc, immediately, on every one of the first 10/10
  //    validation runs. Root cause once reasoned through: this is
  //    structurally the same scheme as attempt 2 (lock v0, v1, then every
  //    ring vertex, unlock-all-and-retry on any single failure), just reached
  //    by a different route -- and address-sorting the acquisition order,
  //    which helped nothing here, only prevents deadlock in a *blocking*
  //    lock-ordering scheme. try_lock_and_get_incident_cells() never blocks,
  //    so ordering has no bearing on the actual failure mode: at internal
  //    flip's candidate density (essentially every internal edge is a
  //    candidate, so ring sets overlap constantly), many candidates
  //    contending for overlapping vertex sets thrash -- fail, unlock, yield,
  //    retry -- and each retry re-pays a fresh Cells_vector allocation per
  //    ring vertex (up to ~8, small_vector<Cell_handle,64> each). That is a
  //    genuine unbounded allocation-churn source under optimistic try-lock at
  //    this contention level, not a bug in the ring-vertex computation
  //    itself. The conclusion attempt 2 already reached stands, now
  //    confirmed via a structurally different code path: optimistic
  //    try-lock-and-retry cannot cheaply protect a multi-vertex footprint
  //    when the candidate graph is this dense, regardless of how the vertex
  //    set is discovered or ordered. A real fix needs either non-optimistic
  //    locking (blocking acquire-in-order, not available on this triangulation's
  //    lock primitive today) or restructuring flip's candidate scheduling so
  //    contending candidates are never scheduled concurrently in the first
  //    place.
  //
  // 8. prefers_boundary_waves(): scope the conflict-free-wave scheme
  //    (attempt 4's mechanism) down to just the BOUNDARY subset instead of
  //    the whole candidate set, so serial wave selection stays proportional
  //    to the (small) boundary set rather than the whole mesh. First cut
  //    (stamping/testing only lv's own 1-ring, like
  //    apply_conflict_free_waves() does for non-halo operations) still
  //    crashed (SIGSEGV) and livelocked (timeouts) on nearly every run:
  //    flip's true footprint needs the same extra cell-layer that
  //    footprint_needs_halo()/Pass 2b exists to cover for classification, and
  //    the wave claim/test wasn't extended to match. Extending both the test
  //    and the claim to that same one-hop-beyond-lv halo (matching Pass 2b's
  //    semantics) fixed the SIGSEGV's specific cause in principle, but the
  //    revalidation run still failed at a similar rate (still one SIGSEGV,
  //    mostly livelock/timeout) -- meaning either a further, undiagnosed gap
  //    remains in the halo-extended claim, or the per-candidate,
  //    per-wave-retry cost of recomputing that extended set is itself now
  //    the bottleneck (each pending candidate re-walks its halo on every wave
  //    pass, not just once). Not fully root-caused before this session ran
  //    out of budget to investigate further -- reverted rather than left in
  //    a partially-diagnosed, still-crashing state. See the "mesh3_243015
  //    flip defect" log entries (ATTEMPT 8 / 8b) for the two validation runs
  //    and the exact rc breakdowns.
  //
  // *** ATTEMPT 8c: re-enabled 8b under crashbt_debug in-process backtrace
  // instrumentation to root-cause the remaining SIGSEGV. Found it:
  // find_best_flip()'s lazy population of the per-worker inc_cells_map cache
  // was racing across threads because the wave's ACCEPTED members were
  // still executed via a nested tbb::parallel_for -- meaning find_best_flip's
  // real read footprint (it explores/caches vertices while evaluating
  // candidate flip rotations) reaches further than any vertex-set claim
  // computed ahead of time could practically capture.
  //
  // *** ATTEMPT 9: execute each wave's accepted members SEQUENTIALLY instead
  // of via tbb::parallel_for (only the deferred boundary subset, not the
  // whole phase). This did eliminate the crash -- 13/13 runs with no
  // SIGSEGV/SIGABRT -- but all 13 timed out (rc=124). CPU usage during a run
  // sat at ~100% (one core), not ~400%, meaning the "boundary" subset this
  // scheme forces onto a single thread is NOT a small minority for internal
  // flip: once locked_vertices() is widened to accurately reflect the real
  // footprint, most internal-edge candidates end up classified boundary
  // (their ring vertices are shared across bucket interfaces often enough,
  // at this candidate density, that "boundary" is closer to "most of the
  // mesh" than "the interface only"). Serializing that subset for
  // correctness therefore serializes most of the flip phase -- not a fixable
  // edge case, a direct consequence of how dense flip's candidate graph is.
  //
  // *** ATTEMPT 9b: before accepting that conclusion, fixed a separate real
  // bug found by inspection: the wave claim's halo-growth step (see
  // Elementary_operation.h) was growing from e.first/e.second's ENTIRE
  // incident-cell star (every cell touching that vertex anywhere in the
  // mesh), not just cells near this specific edge -- a severe over-claim
  // that alone could explain excessive spurious conflicts. Added
  // halo_growth_skip() to skip growing from the edge's own two endpoints,
  // only from its ring vertices (see the base class's comment). Rebuilt,
  // revalidated: still all timeouts, and the live process was still sitting
  // at ~100% CPU -- confirming the over-broad-growth bug was real and worth
  // fixing, but not the (or not the only) cause of the serialization; the
  // "boundary is most of the mesh for this operation" conclusion above
  // stands independent of it.
  //
  // *** CONCLUSION: attempt 8's whole category (claim scoped to a subset
  // rather than the whole mesh) is now genuinely ruled out for internal
  // flip, not just "found buggy in its first implementation" as attempt 7's
  // writeup left it -- the scoping assumption it depends on (boundary is a
  // small minority) does not hold for this operation's candidate density.
  // Reverted (again) to the pre-session baseline: {e.first, e.second} only,
  // no prefers_boundary_waves()/halo_growth_skip() override.
  // *** ATTEMPT 10 (diagnostic): re-enabled with scoped halo_growth_skip()
  // and PARALLEL wave execution (see Elementary_operation.h), to test the
  // one combination not yet tried: scoped growth + parallel, isolating
  // whether the scoping fix alone closes the gap 8b/8c's crash exposed.
  // **It did not** -- run 1/20 crashed (`rc=134`) with the identical
  // find_best_flip/inc_cells_map backtrace signature as 8c's crash (same
  // call chain, offsets shifted only by the intervening code changes). This
  // definitively separates two previously-conflated variables: growth
  // scoping (broad vs. narrow) does not affect whether the crash occurs;
  // only removing PARALLEL wave execution (attempt 9) does. The remaining
  // gap is therefore not a matter of the claim being too narrow or too
  // broad -- it is something a vertex-based claim (however computed) cannot
  // express for this operation. Inspection of `flip_n_to_m`'s cell-adjacency
  // stitching (~line 926-969: `ch->set_neighbor(v, facet.first)` writes into
  // `facet.first`, a MIRRORED neighbor cell one hop across the ring boundary
  // via `tr.mirror_facet()`) is the most likely site: that neighbor cell's
  // own 4th vertex is not necessarily anything the claim stamps directly,
  // even though the neighbor cell itself IS reachable through a ring
  // vertex's incident-cell walk in principle -- suggesting either a subtle
  // gap in that specific reachability argument, or unsynchronized
  // `Concurrent_compact_container`-level cell creation/deletion contention
  // that no per-vertex claim could ever cover. Not confirmed further within
  // this session's remaining budget. Reverted (again) to the pre-session
  // baseline: {e.first, e.second} only, no prefers_boundary_waves()/
  // halo_growth_skip() override.
  // *** ATTEMPT 11 (diagnostic, debug-symbol build): reproduced attempt 10's
  // crash under -g, getting a line-level backtrace instead of function-name-
  // only. New finding: the crash resolves to
  // `Time_stamper<...>::hash_value()` -> `vertex->time_stamp()`
  // (SMDS_3/include/CGAL/Simplicial_mesh_vertex_base_3.h:172), called from
  // `CC_iterator::operator->()` (STL_Extension/include/CGAL/Compact_container.h:1121)
  // -- i.e. the crash is inside HASHING a Vertex_handle key while inserting
  // into the per-worker `inc_cells` unordered_map, not inside the
  // incident-cells walk or the adjacency-stitching write theorized earlier.
  // This points at a genuinely different mechanism than the "unclaimed
  // mirrored-neighbor-cell write" hypothesis from attempt 10: either (a) the
  // Vertex_handle being hashed is itself dangling/stale, or (b) time_stamp()
  // unsafe to read concurrently. CORRECTED (external review, 2026-08-21):
  // flip never deletes vertices (only cells -- tds().delete_cell() at
  // :329/:979, tds().create_cell() at :899), so (b) is dead and (a) was
  // mis-localized. The real mechanism: populating inc_cells[w] for an
  // UNCLAIMED ring vertex w calls tr.incident_cells(vh, ...) (:607, :785,
  // :1789, :2405), which walks w's cell star via a circulator following
  // w->cell() -- unprotected against a *concurrent* set_cell() (:318, :964),
  // delete_cell() (:329, :979), or create_cell() (:899) on exactly those
  // cells, from another worker whose own ring overlaps this one's. The
  // circulator can follow a pointer into freed/recycled
  // Concurrent_compact_container storage, yield a garbage Cell_handle, and
  // ch->vertex(v) off it hands a garbage Vertex_handle straight into the
  // inc_cells[...] insert -- landing on the reported crash frame
  // (hash_value -> CC_iterator::operator-> -> time_stamp()). The crash site
  // is real but downstream; the corruption is in the unprotected star walk.
  //
  // *** ATTEMPT 12/12b: re-attempt #5 (lock_zone_precomputed(), locking the
  // Pass-1-safely-discovered ring-vertex set instead of rediscovering it
  // live), with contention mitigation for its earlier bad_alloc failure
  // (bounded retry + backoff + serial fallback on cap -- 64 attempts in
  // 12, 4 in 12b, see Elementary_operation.h's Pass 3 boundary branch).
  // **Both still hit std::bad_alloc**, identically, on the very first run
  // (`terminate called after throwing an instance of 'std::bad_alloc'`,
  // right after the collapse phase). Tightening the retry cap from 64 to 4
  // made no difference -- this rules out "retry-churn count" as the driver
  // of the allocation pressure (attempt 2's original diagnosis). The cost
  // is per-*attempt*, not cumulative-across-retries: at this candidate
  // density (most candidates are boundary, per attempt 9), many threads
  // simultaneously call try_lock_and_get_incident_cells() for overlapping
  // ring vertices, and if any of those vertices has high valence, even a
  // single attempt's Cells_vector can be large -- multiplied across many
  // concurrently-contending threads, this plausibly reaches genuine
  // multi-GB territory regardless of how few retries each candidate is
  // allowed. The externally-reviewed mechanism (unprotected
  // incident_cells(w) walk racing set_cell()/delete_cell()/create_cell() on
  // an unclaimed ring vertex's cell star) is very likely still the correct
  // explanation for the ORIGINAL crash; what this rules out is "lock every
  // ring vertex, retry on contention" as a *viable implementation* of the
  // fix it implies -- the locking itself is too expensive at this density
  // to survive contention, independent of retry tuning. A real
  // implementation of this mechanism's fix needs to reduce the PER-ATTEMPT
  // cost (e.g. cache/reuse each vertex's incident-cell list across retries
  // instead of re-fetching from scratch every attempt, or find a way to
  // avoid materializing the incident-cell list at lock time at all -- only
  // the *protection* is needed at this stage, not the list itself, which
  // find_best_flip re-derives when it actually needs it). Reverted (again)
  // to the pre-session baseline: {e.first, e.second} only, no
  // prefers_boundary_waves()/halo_growth_skip() override, no
  // lock_zone_precomputed() override (falls back to lock_zone() below).
  // *** ATTEMPT 13: same widened-ring-vertex idea as 12/12b, but
  // lock_zone_precomputed() now only *locks* each vertex (try_lock_vertex(),
  // no cell-list fetch) instead of try_lock_and_get_incident_cells(). The
  // incident-cell cache is left to execute_operation()'s existing lazy
  // population (`if (o_inc_vh.empty()) tr.incident_cells(vh, ...)`,
  // flip_edges.h ~605-607/783-785) -- now safe to run there because it only
  // ever runs after every vertex in `lv` is already locked, so no other
  // locked-path candidate can be concurrently mutating (set_cell/
  // delete_cell/create_cell) any cell incident to those vertices. This
  // directly targets attempt 12/12b's finding: the bad_alloc came from
  // eagerly building a Cells_vector for every vertex on every lock attempt,
  // not from the number of retries -- try_lock_vertex() does no allocation
  // at all, so a failed attempt costs a lock CAS per vertex, not a cell walk.
  // Locks are acquired on demand inside execute_operation(), so they must
  // survive from one element to the next for the incident-cell cache built
  // under them to stay usable. See Elementary_operation::prefers_bucket_scoped_locks().
  bool prefers_bucket_scoped_locks() const override { return true; }

  // The worker just dropped its locks: every cached cell star was derived
  // from the triangulation while those locks were held and may now be stale.
  void on_locks_released() override { inc_cells_map().clear(); }

  void locked_vertices(const Element_type& e, const C3t3& c3t3,
                       boost::container::small_vector<Vertex_handle, 2>& out) const override
  {
    out.push_back(e.first);
    out.push_back(e.second);

    auto& tr = c3t3.triangulation();
    Cell_handle ch; int i0 = 0, i1 = 0;
    if (!tr.tds().is_edge(e.first, e.second, ch, i0, i1))
      return;

    auto circ = tr.incident_cells(Edge(ch, i0, i1));
    const auto done = circ;
    do
    {
      for (int k = 0; k < 4; ++k)
      {
        const Vertex_handle w = circ->vertex(k);
        if (w != e.first && w != e.second
            && std::find(out.begin(), out.end(), w) == out.end())
          out.push_back(w);
      }
    } while (++circ != done);
  }

  // Locks every vertex in `lv`, sorted by address, using the lock-only
  // primitive (no incident-cell fetch -- see the comment above
  // locked_vertices() for why that matters). On any failure, returns false
  // without unlocking; the caller (Elementary_operation.h's Pass 3 boundary
  // branch) unlocks everything acquired so far before retrying.
  bool lock_zone_precomputed(const Element_type&, const C3t3& c3t3,
                              const boost::container::small_vector<Vertex_handle, 2>& lv) const override
  {
    auto& tr = c3t3.triangulation();
    boost::container::small_vector<Vertex_handle, 8> ordered(lv.begin(), lv.end());
    std::sort(ordered.begin(), ordered.end(),
              [](const Vertex_handle& a, const Vertex_handle& b) { return &*a < &*b; });

    for (const Vertex_handle& v : ordered)
      if (!tr.try_lock_vertex(v))
        return false;
    return true;
  }

  bool footprint_needs_halo() const override { return true; }

  bool requires_ordered_processing() const override { return false; }

#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
  // See flip_buckets_per_thread() above : flip holds bucket-scoped locks, so it
  // wants many short-lived buckets where smooth wants few large ones.
  std::size_t buckets_per_thread_hint() const override
  { return flip_buckets_per_thread(); }
#endif

  std::string operation_name() const override { return "Flip edges (internal)"; }
};

// Flip of boundary edges. Mirrors the former flipBoundaryEdges(): compute the
// per-vertex boundary valences once, then, for each boundary edge, flip on the
// surface when it lowers the valence cost, updating the valences accordingly.
template <typename C3t3, typename CellSelector, typename Visitor>
class Boundary_edge_flip_operation
    : public Edge_flip_operation_base<C3t3, CellSelector, Visitor>,
      public Elementary_operation<C3t3,
                                 std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                                 std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>>
{
  using BaseClass = Edge_flip_operation_base<C3t3, CellSelector, Visitor>;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::Facet;
  using typename BaseClass::Cells_vector;
  using BaseClass::m_cell_selector;
  using BaseClass::m_visitor;
  using BaseClass::inc_cells;
  // Dependent base member: GCC/Clang need the using-declaration, MSVC does not.
  using BaseClass::inc_cells_map;

  using Subdomain_index = typename C3t3::Subdomain_index;
  using Surface_patch_index = typename C3t3::Surface_patch_index;
  using Spi_map = boost::unordered_map<Surface_patch_index, unsigned int>;

  mutable boost::unordered_map<Vertex_handle, Spi_map> m_boundary_vertices_valences;
  // execute_operation() runs on Pass 3's parallel_for, where different
  // buckets execute truly concurrently on different worker threads -- unlike
  // the triangulation itself, m_boundary_vertices_valences is a plain,
  // non-thread-safe boost::unordered_map, and flip_surface_edge() both reads
  // it (to decide whether to flip) and writes it (to record the result) in
  // one read-decide-write sequence. That sequence needs to be atomic as a
  // whole -- a per-entry atomic wouldn't preserve the invariant the cost
  // comparison relies on -- so it is serialized here rather than made
  // lock-free. Boundary edges are a minority of flip candidates, so this
  // does not serialize the bulk of the phase, only their processing among
  // themselves. lock_zone()'s triangulation locking does not protect this
  // map at all, so this hazard exists independently of any elision
  // classification (locked_vertices()/footprint_needs_halo() above).
  mutable std::mutex m_valence_mutex;

public:
  using Incident_cells_map = typename BaseClass::Incident_cells_map;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using Base_operation = Elementary_operation<C3t3, Edge_vv, std::vector<Edge_vv>>;
  using Element_type = typename Base_operation::Element_type;
  static_assert(std::is_same_v<Element_type, Edge_vv>, "Element_type must be Edge_vv");
  using Element_range = typename Base_operation::Element_range;

  Boundary_edge_flip_operation(CellSelector& cell_selector, Visitor& visitor, Incident_cells_map& incident_cells)
      : BaseClass(cell_selector, visitor, incident_cells) {}

  Element_range get_elements(const C3t3& c3t3) const override
  {
    std::vector<Edge> boundary_edges;
    boost::unordered_map<Vertex_handle, std::unordered_set<Subdomain_index>> vertices_subdomain_indices;
    m_boundary_vertices_valences.clear();
    collectBoundaryEdgesAndComputeVerticesValences(c3t3,
                                                   m_cell_selector,
                                                   boundary_edges,
                                                   m_boundary_vertices_valences,
                                                   vertices_subdomain_indices);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    if (!debug::are_cell_orientations_valid(c3t3.triangulation()))
      std::cerr << "ERROR in ORIENTATION" << std::endl;
#endif

    std::vector<Edge_vv> candidate_edges_for_flip;
    for (const Edge& e : boundary_edges)
    {
      if (!c3t3.is_in_complex(e))
        candidate_edges_for_flip.push_back(make_vertex_pair(e));
    }
    return candidate_edges_for_flip;
  }

  bool execute_operation(const Element_type& vp, C3t3& c3t3) override
  {
    const Vertex_handle vh0 = vp.first;
    const Vertex_handle vh1 = vp.second;
    typename C3t3::Triangulation& tr = c3t3.triangulation();

    auto& inc_cells_l = inc_cells_map(); // per-worker cache
    Cells_vector& inc_vh0 = inc_cells_l[vh0];
    if (inc_vh0.empty())
      tr.incident_cells(vh0, std::back_inserter(inc_vh0));

    Cell_handle c;
    int i, j;
    if (!is_edge_uv(vh0, vh1, inc_vh0, c, i, j))
      return false;

    Edge edge(c, i, j);
    std::vector<Facet> boundary_facets;
    const bool on_boundary = is_boundary_edge(edge, c3t3, m_cell_selector, boundary_facets);

    if (!on_boundary || boundary_facets.size() != 2)
      return false;

    // flip_surface_edge() reads and then writes m_boundary_vertices_valences
    // as one sequence; see the comment on m_valence_mutex above for why that
    // needs to be atomic as a whole, independent of triangulation locking.
    std::lock_guard<std::mutex> valence_lock(m_valence_mutex);
    return flip_surface_edge(c3t3, edge,
                              boundary_facets,
                              m_boundary_vertices_valences,
                              inc_cells_l,
                              MIN_ANGLE_BASED,
                              m_visitor);
  }

  bool lock_zone(const Element_type& e, const C3t3& c3t3) const override
  {
    auto& tr = c3t3.triangulation();
    Cells_vector inc_cells_first, inc_cells_second;
    const bool locked = tr.try_lock_and_get_incident_cells(e.first, inc_cells_first)
                      && tr.try_lock_and_get_incident_cells(e.second, inc_cells_second);
    if (locked)
    {
      auto& inc_cells_l = inc_cells_map(); // per-worker cache
      inc_cells_l[e.first] = inc_cells_first;
      inc_cells_l[e.second] = inc_cells_second;
    }
    return locked;
  }

  // Unlike Internal_edge_flip_operation, boundary flip's own triangulation
  // footprint never reaches past its edge's two endpoints: flip_on_surface()
  // and flip_n_to_m_on_surface() (this file) mutate only cells_around_edge --
  // cells incident to the edge itself, already inside {e.first, e.second}'s
  // own claim. No ring-vertex widening or halo growth needed. (This method
  // briefly carried a ring-vertex widening copied from
  // Internal_edge_flip_operation by analogy, without checking that
  // assumption against boundary flip's own code -- it doesn't have
  // find_best_flip()'s full-star read. Reverted; the real remaining hazard
  // for boundary flip turned out to be m_boundary_vertices_valences, fixed
  // above with m_valence_mutex, not a cell/vertex-star gap.)
  void locked_vertices(const Element_type& e, const C3t3&,
                       boost::container::small_vector<Vertex_handle, 2>& out) const override
  {
    out = { e.first, e.second };
  }

  bool footprint_needs_halo() const override { return false; }

  bool requires_ordered_processing() const override { return false; }

  std::string operation_name() const override { return "Flip edges (boundary)"; }
};

}//namespace internal
}//namespace Tetrahedral_remeshing
}//namespace CGAL

#endif // CGAL_INTERNAL_FLIP_EDGES_H
