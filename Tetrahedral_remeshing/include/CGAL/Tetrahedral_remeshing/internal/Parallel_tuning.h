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

#include <cstdio>
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
  // A1/A2 -- SHIPPED 2026-09-05. A cell adjacent to a locked star shares a
  // facet with it, so three of its four vertices are already held; only the
  // apex needs locking. Measured independently: +3.33% (collapse) and +6.59%
  // (flip), both 24/24 configs faster, cycles improving at least as much as
  // wall. Verified with the lock probe against a liveness reference of 118,419
  // unlocked writes, and by a 60-run crash soak.
  // Default ON. `CGAL_TR_APEX_HALO_*=0` restores the pre-ship behaviour.
  bool apex_only_collapse_halo = true;
  bool apex_only_flip_halo     = true;
  // A3/A4: hoist the thread-local lock-grid handle out of try_lock, so the
  // enumerable_thread_specific lookup is paid once per zone instead of once
  // per vertex of every cell in it. Two INDEPENDENT regions, so that each can
  // be measured on its own: the halo (A3) and the star walk (A4).
  // MVLZ-SPLIT. Measured 2026-09-05 (MVLZ_SPLIT.md): a split touches only the
  // cells incident to its edge and their facet-neighbours -- 92.5M traced
  // accesses over 10 splits, 0 loads and 0 stores on any other cell, and
  // nothing beyond graph distance 1. lock_zone() nevertheless takes both full
  // endpoint stars, which is 2.1x the vertices. This switch takes the measured
  // set instead: the edge's cell ring, plus one apex per ring facet.
  // Default OFF. This is a correctness-critical change, so it ships only after
  // a crash soak, and the OFF arm is the control in that soak.
  // 0 = today's zone (both full endpoint stars)
  // 1 = the measured MVLZ (edge ring + one apex per ring facet)
  // 2 = SABOTAGE, endpoints only. Deliberately insufficient: the positive
  //     control for the lock-coverage check. A "0 violations" result from
  //     mode 1 means nothing unless the same check is shown to FIRE here.
  int  mvlz_split_zone         = 0;
  // The same three arms for collapse (MVLZ_COLLAPSE.md).
  // 0 = today's zone (both full stars + the A1 apex halo + the midpoint)
  // 1 = the measured MVLZ
  // 2 = SABOTAGE, endpoints + midpoint only. Deliberately insufficient.
  int  mvlz_collapse_zone      = 0;
  // The same arms for flip. NOTE the prior here is the OPPOSITE of collapse's:
  // `zone_ring` is 2 for flip, the class comment says a flip re-stitches mirror
  // cells that "live in the two-ring", and `CGAL_TR_FLIP_HALO_LOCK=0` is
  // documented as crashing. So mode 1 is expected to FAIL, and that expectation
  // is what makes it worth running: a control that is supposed to fire.
  // 0 = today's zone (both full stars + the apex halo)
  // 1 = both stars, no halo -- a CONTROL that is expected to fail, not the
  //     candidate. Measured: it takes unprotected stores from 110 to 310 on
  //     the same 8 operations (MVLZ_FLIP.md 4b)
  // 3 = the measured MVLZ: ring cells + their mirror cells. 13.00 zone
  //     vertices against 32.00 today, a strict subset in 24/24 probed
  //     operations. THE CANDIDATE
  // (2 = sabotage is NOT implemented yet; the loader rejects it rather than
  //  silently behaving like 0 -- an arm that does not differ from its control
  //  is the A-vs-A' failure this campaign has already paid for twelve times.)
  int  mvlz_flip_zone          = 0;
  // N9: give the collapse re-queue a patch cache. Scoped to one operation.
  bool requeue_patch_cache     = false;
  // Private-marking star walk in surface_patch_index(): no shared tds_data.
  bool private_marking         = false;
  // Flip's equivalent of `private_marking`, and NOT covered by it: flip's
  // shared-byte marking comes from `incident_cells()` and `tds().is_edge()`
  // walks started at a ring apex, not from `surface_patch_index()`. Measured
  // 2026-09-07: without this, a flip stores at depth 2 on every mesh and class
  // traced, 72-110 of those stores unprotected WITH the shipped zone in place.
  // Parallel path only. See MVLZ_FLIP.md §4.
  bool private_flip_marking    = false;
  bool halo_tls_hoist          = false;
  bool star_tls_hoist          = false;
  // A5: size the lock grid from mesh density instead of a constant.
  // 0 disables; otherwise the target number of star-sized neighbourhoods per
  // grid cell.
  int  lock_grid_per_star      = 0;
  // A6 -- SHIPPED 2026-09-05, +0.62%. Bounded attempts, then defer to the end
  // of the bucket, instead of spinning on a contended zone. Unordered
  // operations only. The gain is NOT reclaimed spin -- utilisation is
  // unchanged, 3.416 cores against 3.413 -- it is not discarding the two-star
  // BFS walk on each failed attempt. Default ON.
  bool defer_on_conflict       = true;

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
  // D4 -- SHIPPED 2026-09-05, +0.44%. 0 submission order, 1 largest-bucket-
  // first (LPT). Equal-count buckets hold equal counts but not equal work, so
  // the phase waits on its slowest bucket; submitting the big ones first
  // shortens the tail for the cost of one sort of ~16 vectors. Default 1.
  int  bucket_schedule         = 1;

  /**
  * Turns the four changes shipped on 2026-09-05 off together, so their
  * COMBINED effect can be measured against the state that preceded them in a
  * single A/B with one environment variable (POLICY 0.2). `CGAL_TR_SHIP4=0`
  * is the pre-ship arm; unset or 1 is what ships.
  */
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
  // For a switch that is ON by default: unset means on, and an explicit 0
  // turns it off. `flag()` cannot express that.
  static bool flag_on(const char* name)
  {
    const char* const e = std::getenv(name);
    return (e == nullptr) || (std::atoi(e) != 0);
  }
  // For an integer whose 0 is a meaningful value rather than "use the default".
  static int number_or(const char* name, const int dflt)
  {
    const char* const e = std::getenv(name);
    return (e == nullptr) ? dflt : std::atoi(e);
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
    t.apex_only_collapse_halo  = flag_on("CGAL_TR_APEX_HALO_COLLAPSE");
    t.apex_only_flip_halo      = flag_on("CGAL_TR_APEX_HALO_FLIP");
    t.mvlz_split_zone          = number("CGAL_TR_MVLZ_SPLIT_ZONE", 0);
    t.mvlz_collapse_zone       = number("CGAL_TR_MVLZ_COLLAPSE_ZONE", 0);
    t.mvlz_flip_zone           = number("CGAL_TR_MVLZ_FLIP_ZONE", 0);
    if (t.mvlz_flip_zone != 0 && t.mvlz_flip_zone != 1 && t.mvlz_flip_zone != 3)
    {
      std::fprintf(stderr, "CGAL_TR_MVLZ_FLIP_ZONE=%d is not implemented\n",
                   t.mvlz_flip_zone);
      std::abort();
    }
    t.requeue_patch_cache      = flag("CGAL_TR_REQUEUE_PATCH_CACHE");
    t.private_marking          = flag("CGAL_TR_PRIVATE_MARKING");
    t.private_flip_marking     = flag("CGAL_TR_PRIVATE_FLIP_MARKING");
    t.halo_tls_hoist           = flag("CGAL_TR_HALO_TLS_HOIST");
    t.star_tls_hoist           = flag("CGAL_TR_STAR_TLS_HOIST");
    t.lock_grid_per_star       = number("CGAL_TR_LOCK_GRID_PER_STAR", 0);
    t.defer_on_conflict        = flag_on("CGAL_TR_DEFER_ON_CONFLICT");
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
    t.bucket_schedule          = number_or("CGAL_TR_BUCKET_SCHEDULE", 1);

    // ONE switch for the shippable collapse change, because `ab_alloc.sh` takes
    // a single environment variable as the arm. The two halves are inseparable
    // anyway: the smaller zone is only safe once the marking is private
    // (MVLZ_COLLAPSE.md STATUS), so measuring them apart measures nothing that
    // could ship.
    if (flag("CGAL_TR_COLLAPSE_MVLZ"))
    {
      t.private_marking    = true;
      t.mvlz_collapse_zone = 1;
    }

    // ONE switch for the shippable flip change, for the same reason as
    // CGAL_TR_COLLAPSE_MVLZ above: `ab_alloc.sh` takes a single environment
    // variable as the arm, and the two halves are inseparable anyway. The
    // smaller zone is only safe once flip's marking is private -- mode 3 stops
    // locking the stars, and a marking walk over a star it no longer holds is
    // an unprotected write by construction -- so measuring either half alone
    // measures something that cannot ship.
    //
    // NOTE this changed on 2026-09-07: it used to select mode 1 (both stars,
    // no halo), which the measurement says is the WRONG DIRECTION -- it drops
    // protection from the mirror cells, which ARE written, while keeping the
    // rest of the stars, which are not (MVLZ_FLIP.md §5). Mode 1 survives as a
    // control that is supposed to fail.
    if (flag("CGAL_TR_FLIP_MVLZ"))
    {
      t.private_flip_marking = true;
      t.mvlz_flip_zone       = 3;
    }

    /**
    * THE MVLZ ARMS SHIP ON BY DEFAULT from 2026-09-07, in the same shape as
    * CGAL_TR_SHIP4 above: unset or 1 is what ships, `CGAL_TR_MVLZ_SHIP=0` is
    * the pre-ship arm. Both halves of each operation move together because
    * neither is safe alone -- a smaller zone with shared marking is an
    * unprotected write by construction.
    *
    *   collapse  private_marking + mvlz_collapse_zone=1
    *             ACCEPTED 2026-09-07, PGO acceptance, ghat +12.542%,
    *             208 runs, 24/24 configs, Protocol VALID
    *   flip      private_flip_marking + mvlz_flip_zone=3 (ring + mirror)
    *             CLEAR 2026-09-07 on the SCREEN tier, ghat +20.843%,
    *             208 runs, every gate passed, 23/24 configs faster.
    *             POLICY 5.2 makes CLEAR ship-eligible without a PGO run.
    *
    * !! THE A/B VARIABLE CHANGED. !! Before this, `ab_alloc.sh
    * CGAL_TR_COLLAPSE_MVLZ ... 0 1` and `... CGAL_TR_FLIP_MVLZ ... 0 1` were
    * the arms. They are now BOTH ON in either arm, so those A/Bs would compare
    * A against A' and report a null -- the failure this campaign has already
    * paid for twelve times. Use `CGAL_TR_MVLZ_SHIP` as the arm variable, and
    * confirm with envspy that the arms differ before believing any result.
    */
    if (flag_on("CGAL_TR_MVLZ_SHIP"))
    {
      t.private_marking      = true;
      t.mvlz_collapse_zone   = 1;
      t.private_flip_marking = true;
      t.mvlz_flip_zone       = 3;
      // Split joins the shipped set at the user's direction, 2026-09-07.
      // WEAKER EVIDENCE THAN THE OTHER TWO, AND THE DIFFERENCE MATTERS:
      // collapse has a PGO acceptance and flip a screen CLEAR; split has
      // neither. Output is byte-identical to split_zone=0 at one thread on
      // 118287, 124534, 65619 and 102041 f=0.5 (mesh md5), which is a real
      // gate and drift-immune. Its TIMING is not established -- the runs
      // measured ~+1%, inside noise, on a machine that was drifting 25% over
      // the sweep. And its lock-coverage positive control (the `=2` sabotage
      // arm, which must be shown to FIRE before a "0 violations" result from
      // arm 1 means anything) has never been run. Revert this hunk alone if a
      // parallel correctness question appears.
      t.mvlz_split_zone      = 1;
    }

    // One switch for the combined re-measurement.
    if (!flag_on("CGAL_TR_SHIP4"))
    {
      t.apex_only_collapse_halo = false;
      t.apex_only_flip_halo     = false;
      t.defer_on_conflict       = false;
      t.bucket_schedule         = 0;
    }
    return t;
  }
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_PARALLEL_TUNING_H
