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

#ifndef CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_3_H
#define CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_3_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Simplicial_mesh_vertex_base_3.h>

namespace CGAL {
namespace Tetrahedral_remeshing {

/*!
\ingroup PkgTetrahedralRemeshingClasses

The class `Remeshing_vertex_base_3` is a model of the concept `RemeshingVertexBase_3`.
It is designed to serve as vertex base class for the 3D triangulation
used in the tetrahedral remeshing process.

\tparam Gt is the geometric traits class.
It must be a model of the concept `RemeshingTriangulationTraits_3`.

\tparam Vb is a vertex base class from which `Remeshing_vertex_base_3` derives.
It must be a model of the concept `SimplicialMeshVertexBase_3`.

\cgalModels{RemeshingVertexBase_3,SimplicialMeshVertexBase_3}
*/
template<typename Gt,
         typename Vb = CGAL::Simplicial_mesh_vertex_base_3<Gt,
                         int /*Subdomain_index*/,
                         int /*Surface_patch_index*/,
                         int /*Curve_index*/,
                         int /*Corner_index*/> >
class Remeshing_vertex_base_3
  : public Vb
{
public:
  template <typename TDS2>
  struct Rebind_TDS
  {
    using Vb2 = typename Vb::template Rebind_TDS<TDS2>::Other;
    using Other = Remeshing_vertex_base_3<Gt, Vb2>;
  };

public:
  using Vb::Vb; // constructors

  // Transient dense index assigned per smooth phase by the vertex smoother
  // (Vertex_smoothing_context::reset_vertex_id_map). Stored on the vertex to
  // avoid an unordered_map<Vertex_handle, std::size_t> whose CC_iterator hashing
  // and per-lookup cost dominated the serial refresh() bookkeeping. Not
  // serialized; reassigned 0..N-1 over finite_vertex_handles() at the start of
  // every smooth phase, so values from prior iterations are always overwritten.
  std::size_t smoothing_id() const { return smoothing_id_; }
  void set_smoothing_id(const std::size_t id) { smoothing_id_ = id; }

  // Transient stamp marking the vertex as claimed by the wave a parallel
  // executor is currently assembling (see next_wave_claim_stamp()). Stored on
  // the vertex for the same reason as smoothing_id_ : the alternative, a
  // unordered_set<Vertex_handle>, hashes CC_iterators, and selection touches
  // every vertex of every incident cell of every candidate, so that hashing
  // dominated the phase. Stamps come from a single monotonically increasing
  // counter shared by every executor that claims vertices, and are never
  // reset, so a stale value (including one on a recycled vertex) is always
  // below the current wave's and reads as unclaimed. Not serialized.
  std::size_t wave_claim_stamp() const { return wave_claim_stamp_; }
  void set_wave_claim_stamp(const std::size_t s) { wave_claim_stamp_ = s; }

private:
  std::size_t smoothing_id_ = 0;
  std::size_t wave_claim_stamp_ = 0;
};

} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif //CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_3_H
