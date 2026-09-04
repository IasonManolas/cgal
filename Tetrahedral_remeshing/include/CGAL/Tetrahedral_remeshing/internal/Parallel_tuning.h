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

#ifndef CGAL_TETRAHEDRAL_REMESHING_PARALLEL_TUNING_H
#define CGAL_TETRAHEDRAL_REMESHING_PARALLEL_TUNING_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <cstdlib>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

/**
* One registry for every experimental parallel-tuning switch, read from the
* environment ONCE.
*
* Why a registry rather than a static-local `getenv` per switch, which is what
* the earlier candidates used:
*
*  - POLICY 0.2 requires both arms of an A/B to live in ONE binary, selected at
*    run time. With twenty candidates in flight that is twenty switches, and a
*    per-switch function-local static costs a guard variable load on whatever
*    path reads it. Here the whole struct is filled once and hot code holds a
*    `const Parallel_tuning&`, so a switch is a plain load from a hot cache
*    line.
*  - Every candidate ADDS a field. Two candidates therefore never edit the same
*    line, in any order, which is what lets a patch survive the baseline moving
*    under it.
*
* Defaults are the shipped behaviour. Every switch off must reproduce the
* pre-candidate code path exactly.
*/
struct Parallel_tuning
{
  // ---- Group A : lock area ------------------------------------------------
  // A1/A2: a cell adjacent to a locked star shares a facet with it, so three
  // of its four vertices are already held; only the apex needs locking.
  bool apex_only_collapse_halo = false;
  bool apex_only_flip_halo     = false;
  // A3/A4: hoist the thread-local lock-grid handle out of try_lock, so the
  // enumerable_thread_specific lookup is paid once per zone instead of once
  // per vertex of every cell in it. Two INDEPENDENT regions, so that each can
  // be measured on its own: the halo (A3) and the star walk (A4).
  bool halo_tls_hoist          = false;
  bool star_tls_hoist          = false;
  // A5: size the lock grid from mesh density instead of a constant.
  // 0 disables; otherwise the target number of star-sized neighbourhoods per
  // grid cell.
  int  lock_grid_per_star      = 0;
  // A6: bounded attempts, then defer to the end of the bucket, instead of
  // spinning on a contended zone. Unordered operations only.
  bool defer_on_conflict       = false;

  // ---- Group B : lock elision --------------------------------------------
  // 0 off, 1 R19 as shipped, 2 B1 candidate-local, 3 B2 grid ownership,
  // 4 B3 bucket-scoped acquisition, 5 B4 eight-colour sweep.
  int  elision_mode            = 0;

  // ---- Group C : serial fraction -----------------------------------------
  bool fused_edge_pass         = false;
  bool parallel_refresh        = false;
  bool parallel_normals        = false;
  bool parallel_surface_indices= false;
  bool parallel_boundary_flip  = false;
  bool parallel_candidate_sort = false;

  // ---- Group D : partitioning --------------------------------------------
  // 0 equal-count kd (shipped), 1 Hilbert chunks, 2 METIS.
  int  partitioner             = 0;
  bool reuse_partition         = false;
  int  buckets_per_thread      = 4;
  // 0 submission order, 1 largest-bucket-first (LPT).
  int  bucket_schedule         = 0;

  static const Parallel_tuning& get()
  {
    static const Parallel_tuning t = load();
    return t;
  }

private:
  static bool flag(const char* name)
  {
    const char* const e = std::getenv(name);
    return (e != nullptr) && (std::atoi(e) != 0);
  }
  static int number(const char* name, const int dflt)
  {
    const char* const e = std::getenv(name);
    if (e == nullptr)
      return dflt;
    const int v = std::atoi(e);
    return (v > 0) ? v : dflt;
  }

  static Parallel_tuning load()
  {
    Parallel_tuning t;
    t.apex_only_collapse_halo  = flag("CGAL_TR_APEX_HALO_COLLAPSE");
    t.apex_only_flip_halo      = flag("CGAL_TR_APEX_HALO_FLIP");
    t.halo_tls_hoist           = flag("CGAL_TR_HALO_TLS_HOIST");
    t.star_tls_hoist           = flag("CGAL_TR_STAR_TLS_HOIST");
    t.lock_grid_per_star       = number("CGAL_TR_LOCK_GRID_PER_STAR", 0);
    t.defer_on_conflict        = flag("CGAL_TR_DEFER_ON_CONFLICT");
    t.elision_mode             = number("CGAL_TR_ELISION_MODE", 0);
    t.fused_edge_pass          = flag("CGAL_TR_FUSED_EDGE_PASS");
    t.parallel_refresh         = flag("CGAL_TR_PARALLEL_REFRESH");
    t.parallel_normals         = flag("CGAL_TR_PARALLEL_NORMALS");
    t.parallel_surface_indices = flag("CGAL_TR_PARALLEL_SURF_IDX");
    t.parallel_boundary_flip   = flag("CGAL_TR_PARALLEL_BFLIP_SCAN");
    t.parallel_candidate_sort  = flag("CGAL_TR_PARALLEL_SORT");
    t.partitioner              = number("CGAL_TR_PARTITIONER", 0);
    t.reuse_partition          = flag("CGAL_TR_REUSE_PARTITION");
    t.buckets_per_thread       = number("CGAL_TR_BUCKETS_PER_THREAD", 4);
    t.bucket_schedule          = number("CGAL_TR_BUCKET_SCHEDULE", 0);
    return t;
  }
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_PARALLEL_TUNING_H
