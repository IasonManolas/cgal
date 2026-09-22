// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#ifndef CGAL_INTERNAL_TETRAHEDRAL_REMESHING_INSTRUMENTATION_H
#define CGAL_INTERNAL_TETRAHEDRAL_REMESHING_INSTRUMENTATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>

// Counters used to investigate the parallel remeshing, compiled out unless
// their macro is defined. Nothing here is part of the algorithm: it exists so
// that a change to the parallel machinery can be explained, not merely timed.
//
// They are DIAGNOSTICS and never gates. Each one moves with whatever is being
// tested, and none of them predicts wall time on its own -- a zone change that
// cut retries by a quarter still ran slower, because the saving was spent
// elsewhere. Read them next to a time measurement, never instead of one.
//
//   CGAL_TR_LOCKCOUNT   lock-zone retries per elementary operation, i.e. the
//                       spin: how often a worker had to release its zone and
//                       start over because another thread held part of it.
//
//   CGAL_TR_TOPSTAGE    wall and process-CPU time of every top-level stage of
//                       remesh(), so that the run adds up: the six elementary
//                       operations, the serial spatial-sort rebuild inside
//                       split(), the smoothing's refresh(), the edge scan of
//                       resolution_reached(), and postprocess(). A stage whose
//                       wall time does not fall as threads are added is a
//                       serial section, and the sum against the total says how
//                       much of the run no stage accounts for.

#ifdef CGAL_TR_LOCKCOUNT
#include <cstdint>
#include <iostream>
#include <mutex>
#endif

#ifdef CGAL_TR_TOPSTAGE
#include <CGAL/Real_timer.h>
#include <array>
#include <iostream>
#include <map>
#include <string>
#include <ctime>
#endif

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{

/**
* How often `apply_one()` had to give the zone back and try again -- the spin
* a lock zone pays for, counted rather than timed.
*
* Thread-local, because an atomic in this loop would add exactly the kind of
* contention it is meant to measure; the per-thread blocks are chained into a
* list at first use and summed when the program exits.
*
* A DIAGNOSTIC, not a gate: it moves with the treatment by construction, and a
* crashed run prints nothing at all, so whatever drives it must record the exit
* status separately or a dead run reads as a silent pass.
*/
struct Lockcount_counters
{
  std::uint64_t ops = 0;
  std::uint64_t retries = 0;
  Lockcount_counters* next = nullptr;
  Lockcount_counters();
};

inline std::mutex& lockcount_mutex()
{
  static std::mutex m;
  return m;
}

inline Lockcount_counters*& lockcount_head()
{
  static Lockcount_counters* head = nullptr;
  return head;
}

struct Lockcount_reporter
{
  ~Lockcount_reporter()
  {
    std::uint64_t ops = 0, retries = 0;
    std::lock_guard<std::mutex> g(lockcount_mutex());
    for (const Lockcount_counters* c = lockcount_head(); c != nullptr; c = c->next)
    {
      ops += c->ops;
      retries += c->retries;
    }
    std::cout << "LOCKCOUNT ops=" << ops << " retries=" << retries
              << " retries_per_op="
              << (ops ? double(retries) / double(ops) : 0.0) << std::endl;
  }
};

inline Lockcount_reporter& lockcount_reporter()
{
  static Lockcount_reporter r;
  return r;
}

inline Lockcount_counters::Lockcount_counters()
{
  std::lock_guard<std::mutex> g(lockcount_mutex());
  lockcount_reporter();          // must outlive every thread-local block
  next = lockcount_head();
  lockcount_head() = this;
}

inline Lockcount_counters& lockcount_counters()
{
  static thread_local Lockcount_counters c;
  return c;
}

#ifdef CGAL_TR_TOPSTAGE
/**
* Wall and process-CPU time of each top-level stage of `remesh()`.
*
* A DIAGNOSTIC, not a gate. It exists to separate a stage that gets faster
* when threads are added from one that does not: the second kind is what caps
* the speedup, and no aggregate wall time can tell the two apart. Process CPU
* is recorded next to wall so that a stage running at one core shows up as
* such rather than being inferred.
*
* Single-threaded accounting on purpose -- every stage is entered by the
* calling thread only, so no lock is needed and the timer adds nothing to the
* parallel region it brackets.
*/
struct Top_stage_times
{
  std::map<std::string, std::array<double, 3>> by_stage;  // wall, cpu, calls
  ~Top_stage_times()
  {
    for (const auto& [name, v] : by_stage)
      std::cout << "TOPSTAGE \"" << name << "\" wall=" << v[0]
                << " cpu=" << v[1] << " calls=" << v[2] << std::endl;
  }
};

inline Top_stage_times& top_stage_times()
{
  static Top_stage_times t;
  return t;
}

inline double top_stage_cpu_seconds()
{
  struct timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  return double(ts.tv_sec) + 1e-9 * double(ts.tv_nsec);
}

struct Top_stage_scope
{
  std::array<double, 3>& slot;
  CGAL::Real_timer t;
  double c0;

  explicit Top_stage_scope(const char* name)
    : slot(top_stage_times().by_stage[name]), c0(top_stage_cpu_seconds())
  {
    t.start();
  }
  ~Top_stage_scope()
  {
    t.stop();
    slot[0] += t.time();
    slot[1] += top_stage_cpu_seconds() - c0;
    slot[2] += 1.0;
  }
};

// The name carries the line number, so that two stages may be bracketed in the
// same block -- remesh() times itself and its finalize() side by side.
#define CGAL_TR_TOPSTAGE_JOIN2(a, b) a##b
#define CGAL_TR_TOPSTAGE_JOIN(a, b) CGAL_TR_TOPSTAGE_JOIN2(a, b)
#define CGAL_TR_TOPSTAGE_SCOPE(NAME) \
  ::CGAL::Tetrahedral_remeshing::internal::Top_stage_scope \
    CGAL_TR_TOPSTAGE_JOIN(cgal_tr_topstage_scope_, __LINE__)((NAME))
#else
#define CGAL_TR_TOPSTAGE_SCOPE(NAME) do {} while (0)
#endif

} // end namespace internal
} // end namespace Tetrahedral_remeshing
} // end namespace CGAL

#endif // CGAL_INTERNAL_TETRAHEDRAL_REMESHING_INSTRUMENTATION_H
