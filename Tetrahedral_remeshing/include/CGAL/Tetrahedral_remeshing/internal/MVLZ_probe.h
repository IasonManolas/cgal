// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// MVLZ probe -- GOAL 1: the minimum viable lock zone of an elementary
// operation, measured rather than argued.
//
//     min_zone(A) = write(A)  u  ( read(A) n W ),   W = U write(B) over all B
//
// ---------------------------------------------------------------------------
// TWO INSTRUMENTS, ON PURPOSE
// ---------------------------------------------------------------------------
// TRACE mode (CGAL_TR_MVLZ_TRACE=<manifest path>) is the real instrument. It
// makes the run legible to Valgrind Lackey, which records EVERY load, store
// and modify the operation performs. Nothing here decides what counts as an
// access, so nothing here can miss one -- which matters, because reading the
// source and asserting a footprint is the mode that produced every design
// error in this campaign. It needs two things Lackey cannot supply:
//
//   * SEGMENTATION. Lackey emits addresses, not operation boundaries. This
//     file stores to two magic globals, mvlz_mark_begin / mvlz_mark_end; the
//     post-processor splits the trace on those two addresses. Operations are
//     emitted in program order and the run is single-threaded, so the k-th
//     begin-marker in the trace is the k-th block in the manifest.
//
//   * IDENTITY. A traced address is just a number. Before each probed
//     operation this file writes a manifest block naming every cell and
//     vertex in a ball around the element, with its address, its size, and
//     its BFS depth from the element. A traced address then resolves to
//     object + field offset + graph distance, which is what the question is
//     actually about.
//
// DIFF mode (default) snapshots every byte of every object in the same ball,
// runs the operation, and diffs. It sees only writes, and only NET writes,
// so it cannot answer the question on its own. It is kept because it is
// complete by construction and completely independent of the trace: if the
// diff and Lackey's store set disagree, one of them is wrong and that is
// worth knowing before either is believed.
//
// ---------------------------------------------------------------------------
// HOW DISTANCE IS DEFINED
// ---------------------------------------------------------------------------
// BFS depth on the FINITE vertex adjacency graph from the element's own
// vertices (depth 0). Two rules that decide the result:
//
//   * the infinite vertex is never traversed. It is adjacent to every
//     boundary vertex, so traversing it would put most of the mesh at depth 2
//     and quietly destroy the measurement. Infinite cells stay in the region
//     and are labelled by their finite vertices; accesses to the infinite
//     vertex itself are counted in their own bucket, since no spatial lock
//     can cover an object with no position.
//
//   * a written CELL costs the zone all four of its vertices, because the
//     spatial-lock protocol protects a cell by holding its vertices. So a
//     cell's reported depth is the MAX over its finite vertices, not the min.
//
// ---------------------------------------------------------------------------
// SATURATION
// ---------------------------------------------------------------------------
// The ball has radius CGAL_TR_MVLZ_RADIUS (env, default 2). An access at
// exactly that radius means the ball was too small: a measured maximum equal
// to the search bound is not a measurement, and both modes report INVALID
// rather than a result.
#ifndef CGAL_TETRAHEDRAL_REMESHING_MVLZ_PROBE_H
#define CGAL_TETRAHEDRAL_REMESHING_MVLZ_PROBE_H

#ifdef CGAL_TR_MVLZ_PROBE

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

// ---------------------------------------------------------------------------
// The two magic globals the trace is segmented on. volatile so the stores
// survive optimisation; they are the only reason these variables exist.
// ---------------------------------------------------------------------------
inline volatile unsigned long mvlz_mark_begin = 0;
inline volatile unsigned long mvlz_mark_end   = 0;

struct Mvlz_config
{
  int  radius   = 2;
  long stride   = 1;      // probe every Nth operation
  long max_ops  = -1;     // stop probing after this many (-1 = no limit)
  bool trace    = false;
  bool selftest = false;
  // Emit EVERY cell and vertex of the triangulation in the manifest, not just
  // the ball, tagged with depth OUTSIDE. Then an access that resolves to
  // nothing is PROVABLY not a mesh object, instead of merely unrecognised,
  // and a write beyond the ball shows up as OUTSIDE rather than vanishing.
  // Only affordable on a small mesh, which is the point of using one.
  bool all_objects = false;
  // Stop the process the moment the probed operations are done. Lackey traces
  // the WHOLE run, so the cost is set by how long the program lives, not by
  // how many operations are probed: without this the sphere run spent minutes
  // tracing remeshing we had already finished measuring. _exit() so no
  // destructor work is traced after the last window closes.
  long exit_after = -1;
  // LOCK-COVERAGE CHECK. For every object the operation actually changed,
  // assert this thread held the lock the protocol requires: all four vertices
  // of a changed cell, and the vertex itself for a changed vertex. Combined
  // with the snapshot/diff -- which finds changed objects without knowing
  // where any write happens -- this tests sufficiency of the lock zone
  // COMPLETE BY CONSTRUCTION, instead of at hand-placed write sites. Split
  // has exactly one of the 28 existing sites, so the site-based probe cannot
  // answer this question at all.
  bool lock_check = false;
  // Probe ONE operation kind. Split and collapse are both hooked now, and
  // they share this file's per-instantiation counter, so without a filter
  // CGAL_TR_MVLZ_MAXOPS is a budget split between two kinds by whatever the
  // remesher happened to run first -- cdt/51492 at f=1.5 spent all 400 on
  // splits and reported no collapse at all. A trace with both kinds in it is
  // worse than useless: the post-processor treats every window alike, so the
  // footprint would be the union of two different operations.
  std::string kind_filter;
  // Probe only operations of ONE class (see Mvlz_probe::classify).
  //
  // Needed because sampling the first N operations is BIASED: the collapse
  // work list is shortest-edge-first, the shortest edges sit on surfaces, and
  // 24 of the first 24 collapses on 118287 f=1.5 were non-interior against a
  // whole-run rate of 10.7%. A claim of the form "operations of class X never
  // do Y" cannot be tested on a sample that contains no class X.
  int class_filter = -1;        // -1 = probe every class
  std::FILE* manifest = nullptr;

  static Mvlz_config& get()
  {
    static Mvlz_config c = [] {
      Mvlz_config x;
      if (const char* r = std::getenv("CGAL_TR_MVLZ_RADIUS")) x.radius = std::atoi(r);
      if (const char* s = std::getenv("CGAL_TR_MVLZ_STRIDE")) x.stride = std::atol(s);
      if (const char* m = std::getenv("CGAL_TR_MVLZ_MAXOPS")) x.max_ops = std::atol(m);
      if (const char* t = std::getenv("CGAL_TR_MVLZ_SELFTEST")) x.selftest = (*t == '1');
      if (const char* t = std::getenv("CGAL_TR_MVLZ_ALLOBJ")) x.all_objects = (*t == '1');
      if (const char* t = std::getenv("CGAL_TR_MVLZ_EXIT_AFTER")) x.exit_after = std::atol(t);
      if (const char* t = std::getenv("CGAL_TR_MVLZ_LOCKCHECK")) x.lock_check = (*t == '1');
      if (const char* t = std::getenv("CGAL_TR_MVLZ_KIND")) x.kind_filter = t;
      if (const char* t = std::getenv("CGAL_TR_MVLZ_CLASS")) x.class_filter = std::atoi(t);
      if (const char* p = std::getenv("CGAL_TR_MVLZ_TRACE"))
      {
        x.manifest = std::fopen(p, "w");
        x.trace = (x.manifest != nullptr);
      }
      if (x.stride < 1) x.stride = 1;
      return x;
    }();
    return c;
  }
};

// The trace can be piped straight out of Lackey rather than written to disk,
// and then the post-processor must never block waiting for this file: the
// process writing it is the process whose trace is filling the pipe. So the
// segmentation header is emitted at STATIC INIT, before any mesh work, which
// lets the post-processor open the manifest and start draining the pipe while
// the traced program is still starting up. Emitting it at the first probed
// operation instead deadlocks -- observed: Lackey sat at 0% CPU with an empty
// manifest. The two marker addresses are plain globals and are known this
// early; sizes and radius are template- and config-dependent and follow with
// the first operation block.
inline const bool mvlz_marks_header = []
{
  Mvlz_config& c = Mvlz_config::get();
  if (c.manifest)
  {
    std::fprintf(c.manifest, "#MARKS begin=%p end=%p\n",
                 (void*)&mvlz_mark_begin, (void*)&mvlz_mark_end);
    std::fflush(c.manifest);
  }
  return true;
}();

struct Mvlz_totals
{
  struct Kind
  {
    long ops = 0, ops_written = 0;
    std::vector<long> radius = std::vector<long>(16, 0);
    long saturated = 0;
    long cells_written = 0, cells_created = 0, cells_destroyed = 0;
    long verts_written = 0, verts_created = 0, verts_destroyed = 0;
    long inf_vertex_written = 0;
    long wc_ring = 0, wc_mirror = 0, wc_star_other = 0, wc_far = 0;
    long zone_now_verts = 0, zone_min_verts = 0;
    // lock-coverage results
    long lc_cells_checked = 0, lc_cells_uncovered = 0;
    long lc_verts_checked = 0, lc_verts_uncovered = 0;
    long lc_infinite_skipped = 0, lc_ops_with_violation = 0;
    // split the violations by cause, because the number is useless otherwise
    long lc_uncov_created = 0, lc_uncov_modified = 0;
    long lc_uncov_by_new_vertex = 0, lc_uncov_by_old_vertex = 0;
    long lc_no_lock_ds = 0;
    long f_vertex = 0, f_neighbor = 0, f_point = 0, f_cell = 0;
    std::map<int, long> cell_offsets, vert_offsets;
  };
  std::map<std::string, Kind> kinds;

  static Mvlz_totals& get() { static Mvlz_totals t; return t; }

  void report(const char* path) const
  {
    std::FILE* f = std::fopen(path, "w");
    if (!f) return;
    std::fprintf(f, "# MVLZ diff mode -- NET WRITE footprint by snapshot/diff\n");
    std::fprintf(f, "# radius=%d stride=%ld\n\n",
                 Mvlz_config::get().radius, Mvlz_config::get().stride);
    for (const auto& kv : kinds)
    {
      const Kind& k = kv.second;
      std::fprintf(f, "== %s\nops=%ld ops_that_wrote=%ld\n",
                   kv.first.c_str(), k.ops, k.ops_written);
      if (k.saturated)
        std::fprintf(f, "*** INVALID: %ld ops wrote at the ball bound "
                        "(radius %d); raise CGAL_TR_MVLZ_RADIUS\n",
                     k.saturated, Mvlz_config::get().radius);
      std::fprintf(f, "write-radius histogram (max depth of any net write):\n");
      for (std::size_t i = 0; i < k.radius.size(); ++i)
        if (k.radius[i]) std::fprintf(f, "  r=%zu : %ld\n", i, k.radius[i]);
      std::fprintf(f, "cells written=%ld created=%ld destroyed=%ld\n",
                   k.cells_written, k.cells_created, k.cells_destroyed);
      std::fprintf(f, "verts written=%ld created=%ld destroyed=%ld "
                      "infinite_vertex_written=%ld\n", k.verts_written,
                   k.verts_created, k.verts_destroyed, k.inf_vertex_written);
      std::fprintf(f, "written cells by class: ring=%ld mirror=%ld "
                      "star_other=%ld far=%ld\n", k.wc_ring, k.wc_mirror,
                   k.wc_star_other, k.wc_far);
      if (k.ops)
        std::fprintf(f, "zone vertices per op: locked_today=%.2f minimum=%.2f "
                        "ratio=%.3f\n", double(k.zone_now_verts) / k.ops,
                     double(k.zone_min_verts) / k.ops,
                     k.zone_now_verts ? double(k.zone_min_verts)
                                        / double(k.zone_now_verts) : 0.0);
      std::fprintf(f, "fields moved: vertex()=%ld neighbor()=%ld point()=%ld "
                      "cell()=%ld\n", k.f_vertex, k.f_neighbor, k.f_point, k.f_cell);
      std::fprintf(f, "LOCK COVERAGE of every changed object:\n");
      std::fprintf(f, "  changed cells checked   = %ld   NOT covered = %ld\n",
                   k.lc_cells_checked, k.lc_cells_uncovered);
      std::fprintf(f, "  changed vertices checked= %ld   NOT covered = %ld\n",
                   k.lc_verts_checked, k.lc_verts_uncovered);
      std::fprintf(f, "  infinite vertices skipped (no position, unlockable) = %ld\n",
                   k.lc_infinite_skipped);
      std::fprintf(f, "  operations with >=1 violation = %ld\n", k.lc_ops_with_violation);
      std::fprintf(f, "  uncovered cells: created-by-this-op=%ld  modified-in-place=%ld\n",
                   k.lc_uncov_created, k.lc_uncov_modified);
      std::fprintf(f, "  unheld vertex was: created-by-this-op=%ld  pre-existing=%ld\n",
                   k.lc_uncov_by_new_vertex, k.lc_uncov_by_old_vertex);
      if (k.lc_no_lock_ds)
        std::fprintf(f, "  *** %ld ops had NO lock data structure: the check could not "
                        "run. Build with Parallel_tag.\n", k.lc_no_lock_ds);
      if (k.lc_cells_uncovered || k.lc_verts_uncovered)
        std::fprintf(f, "  *** ZONE INSUFFICIENT: the operation changed an object this "
                        "thread did not hold.\n");
      std::fprintf(f, "changed cell byte offsets:\n");
      for (const auto& o : k.cell_offsets) std::fprintf(f, "  +%d : %ld\n", o.first, o.second);
      std::fprintf(f, "changed vertex byte offsets:\n");
      for (const auto& o : k.vert_offsets) std::fprintf(f, "  +%d : %ld\n", o.first, o.second);
      std::fprintf(f, "\n");
    }
    std::fclose(f);
  }
};

template <typename Tr>
class Mvlz_probe
{
  using Vertex_handle = typename Tr::Vertex_handle;
  using Cell_handle   = typename Tr::Cell_handle;
  using Cell   = typename Tr::Triangulation_data_structure::Cell;
  using Vertex = typename Tr::Triangulation_data_structure::Vertex;

  struct CellSnap
  {
    Cell_handle h; int depth; int klass;
    Vertex_handle v[4]; Cell_handle n[4];
    std::vector<unsigned char> bytes;
  };
  struct VertSnap
  {
    Vertex_handle h; int depth; bool infinite;
    Cell_handle cell; double x, y, z;
    std::vector<unsigned char> bytes;
  };

  Tr& m_tr;
  std::string m_kind;
  bool m_active = false;
  Vertex_handle m_v0, m_v1;
  // The two-zone predicate, recorded per operation so the post-processor can
  // PARTITION the footprint by it instead of pooling every window together.
  // A pooled trace cannot distinguish "this operation class never reaches
  // depth 2" from "the average operation rarely does".
  int m_class = -1;
  std::unordered_map<const Vertex*, int> m_depth;
  std::vector<CellSnap> m_cells;
  std::vector<VertSnap> m_verts;
  std::unordered_set<const Cell*> m_pre_cells, m_ring, m_mirror;
  std::unordered_set<const Vertex*> m_pre_verts;
  long m_zone_now = 0;
  int  m_R = 2;
  // Split's zone is 1-ring(v0) u 1-ring(v1) and nothing else. Collapse's also
  // takes one apex per star facet (the shipped A1 halo). The "locked today"
  // set is what the subset-of-today's-zone validity check compares against, so
  // reporting split's zone for a collapse would make that check compare the
  // measurement against a zone the code does not take.
  bool m_zone_has_halo = false;

  static long& counter() { static long c = 0; return c; }
  static long& windows_closed() { static long c = 0; return c; }
  static bool& header_done() { static bool b = false; return b; }

  int depth_of(Vertex_handle v) const
  {
    auto it = m_depth.find(&*v);
    return (it == m_depth.end()) ? (m_R + 1) : it->second;
  }
  int cell_depth(Cell_handle c) const
  {
    int d = 0;
    for (int i = 0; i < 4; ++i)
      if (!m_tr.is_infinite(c->vertex(i)))
        d = (std::max)(d, depth_of(c->vertex(i)));
    return d;
  }

  // Does this thread hold what the protocol demands for this object?
  // Cell  -> all four of its vertices. Vertex -> itself.
  // The infinite vertex has no position, so no spatial lock can stand for it;
  // it is counted separately rather than called a violation.
  bool cell_covered(Cell_handle c, long& skipped) const
  {
    auto* lds = m_tr.get_lock_data_structure();
    if (lds == nullptr) return true;
    for (int i = 0; i < 4; ++i)
    {
      Vertex_handle v = c->vertex(i);
      if (v == Vertex_handle()) continue;
      if (m_tr.is_infinite(v)) { ++skipped; continue; }
      if (!lds->is_locked_by_this_thread(v->point())) return false;
    }
    return true;
  }
  bool vertex_covered(Vertex_handle v, long& skipped) const
  {
    auto* lds = m_tr.get_lock_data_structure();
    if (lds == nullptr) return true;
    if (v == Vertex_handle()) return true;
    if (m_tr.is_infinite(v)) { ++skipped; return true; }
    return lds->is_locked_by_this_thread(v->point());
  }

  // Lock coverage AT THE MOMENT THE WINDOW OPENS, for the read-side check.
  //   cell   -> 2 = all four vertices held (a writer's requirement)
  //             1 = at least one held      (a reader's requirement: a writer
  //                                         would need all four, so holding
  //                                         one excludes every writer)
  //             0 = none held              (a read here is unprotected)
  //   vertex -> 1 = held. A writer of a vertex holds only that vertex
  //             (smoothing writes point()), so a reader needs it too.
  // The infinite vertex has no position and no spatial lock can stand for it;
  // it counts as held so it does not drown the signal, and is reported apart.
  int cell_coverage(Cell_handle c) const
  {
    auto* lds = m_tr.get_lock_data_structure();
    if (lds == nullptr) return 2;
    int held = 0, n = 0;
    for (int i = 0; i < 4; ++i)
    {
      Vertex_handle v = c->vertex(i);
      if (v == Vertex_handle()) continue;
      ++n;
      if (m_tr.is_infinite(v) || lds->is_locked_by_this_thread(v->point())) ++held;
    }
    if (held == 0) return 0;
    return (held == n) ? 2 : 1;
  }
  int vertex_coverage(Vertex_handle v) const
  {
    auto* lds = m_tr.get_lock_data_structure();
    if (lds == nullptr) return 1;
    if (v == Vertex_handle() || m_tr.is_infinite(v)) return 1;
    return lds->is_locked_by_this_thread(v->point()) ? 1 : 0;
  }

public:
  // The vertex set lock_zone() acquires today, as a set rather than a count,
  // so a smaller measured set can be checked to be a SUBSET of it and not
  // merely smaller.
  void zone_today(std::unordered_set<const Vertex*>& z) const
  {
    std::vector<Cell_handle> inc;
    std::vector<Cell_handle> star;
    for (Vertex_handle v : {m_v0, m_v1})
    {
      inc.clear();
      m_tr.incident_cells(v, std::back_inserter(inc));
      for (Cell_handle c : inc)
      {
        star.push_back(c);
        for (int i = 0; i < 4; ++i) z.insert(&*c->vertex(i));
      }
    }
    if (!m_zone_has_halo) return;
    for (Cell_handle c : star)
      for (int i = 0; i < 4; ++i)
      {
        const Cell_handle n = c->neighbor(i);
        z.insert(&*n->vertex(n->index(c)));      // the apex only -- A1
      }
  }

public:
  explicit Mvlz_probe(Tr& tr) : m_tr(tr) { m_R = Mvlz_config::get().radius; }

  // Call before begin() for an operation whose shipped zone includes the
  // apex halo (collapse), so "locked today" names the zone the code takes.
  void zone_today_has_apex_halo() { m_zone_has_halo = true; }

  // Call before begin(). 1 = the operation is in the class the small zone is
  // claimed to cover, 0 = it is not, -1 = unclassified.
  void classify(int k) { m_class = k; }

  void begin(const char* kind, Vertex_handle a, Vertex_handle b)
  {
    Mvlz_config& cfg = Mvlz_config::get();
    m_kind = kind;
    m_v0 = a; m_v1 = b;
    m_active = false;
    // The filter is applied BEFORE the counter, so stride and max_ops are
    // budgets for the selected kind rather than for whatever ran first.
    if (!cfg.kind_filter.empty() && cfg.kind_filter != m_kind) return;
    if (cfg.class_filter >= 0 && cfg.class_filter != m_class) return;
    const long id = counter()++;
    if (cfg.stride > 1 && (id % cfg.stride) != 0) return;
    if (cfg.max_ops >= 0 && (id / cfg.stride) >= cfg.max_ops) return;
    if (m_tr.is_infinite(a) || m_tr.is_infinite(b)) return;

    m_depth.clear(); m_cells.clear(); m_verts.clear();
    m_pre_cells.clear(); m_pre_verts.clear(); m_ring.clear(); m_mirror.clear();

    // ---- BFS on the finite vertex graph; depth 0 = the element ----
    std::vector<Vertex_handle> frontier{a, b}, next;
    m_depth[&*a] = 0; m_depth[&*b] = 0;
    std::vector<Cell_handle> inc;
    for (int d = 1; d <= m_R + 1 && !frontier.empty(); ++d)
    {
      next.clear();
      for (Vertex_handle v : frontier)
      {
        inc.clear();
        m_tr.incident_cells(v, std::back_inserter(inc));
        for (Cell_handle c : inc)
          for (int i = 0; i < 4; ++i)
          {
            Vertex_handle w = c->vertex(i);
            if (m_tr.is_infinite(w)) continue;        // never traverse infinity
            if (m_depth.emplace(&*w, d).second) next.push_back(w);
          }
      }
      frontier.swap(next);
    }

    // ---- region cells: incident to any vertex at depth <= R ----
    std::vector<Vertex_handle> region_v;
    for (const auto& kv : m_depth)
      if (kv.second <= m_R)
        region_v.push_back(Vertex_handle(const_cast<Vertex*>(kv.first)));

    std::unordered_set<const Cell*> seen;
    for (Vertex_handle v : region_v)
    {
      inc.clear();
      m_tr.incident_cells(v, std::back_inserter(inc));
      for (Cell_handle c : inc)
        if (seen.insert(&*c).second)
        {
          CellSnap s;
          s.h = c; s.depth = cell_depth(c); s.klass = 3;
          for (int i = 0; i < 4; ++i) { s.v[i] = c->vertex(i); s.n[i] = c->neighbor(i); }
          if (!cfg.trace)
            s.bytes.assign(reinterpret_cast<const unsigned char*>(&*c),
                           reinterpret_cast<const unsigned char*>(&*c) + sizeof(Cell));
          m_cells.push_back(std::move(s));
          m_pre_cells.insert(&*c);
        }
    }

    std::unordered_set<const Vertex*> vseen;
    for (const auto& s : m_cells)
      for (int i = 0; i < 4; ++i)
      {
        Vertex_handle w = s.v[i];
        if (!vseen.insert(&*w).second) continue;
        VertSnap t;
        t.h = w;
        t.infinite = m_tr.is_infinite(w);
        t.depth = t.infinite ? -1 : depth_of(w);
        t.cell = w->cell();
        if (!t.infinite) { t.x = w->point().x(); t.y = w->point().y(); t.z = w->point().z(); }
        else             { t.x = t.y = t.z = 0; }
        if (!cfg.trace)
          t.bytes.assign(reinterpret_cast<const unsigned char*>(&*w),
                         reinterpret_cast<const unsigned char*>(&*w) + sizeof(Vertex));
        m_verts.push_back(std::move(t));
        m_pre_verts.insert(&*w);
      }

    // ---- classification: the edge's cell ring, and the ring's mirrors ----
    for (const auto& s : m_cells)
    {
      bool h0 = false, h1 = false;
      for (int i = 0; i < 4; ++i)
      { if (s.v[i] == m_v0) h0 = true; if (s.v[i] == m_v1) h1 = true; }
      if (h0 && h1) m_ring.insert(&*s.h);
    }
    for (const auto& s : m_cells)
      if (m_ring.count(&*s.h))
        for (int i = 0; i < 4; ++i)
          if (!m_ring.count(&*s.n[i])) m_mirror.insert(&*s.n[i]);

    // ---- what is locked today ----
    {
      std::unordered_set<const Vertex*> z;
      zone_today(z);
      m_zone_now = (long)z.size();
    }

    m_active = true;

    // ---- manifest, then the begin marker, then nothing until the op ----
    // Order matters: everything this probe does must land OUTSIDE the marked
    // window, or the probe's own reads would be attributed to the operation.
    if (cfg.trace)
    {
      std::FILE* f = cfg.manifest;
      if (!header_done())
      {
        header_done() = true;
        // #MARKS is already out, written at static init -- see mvlz_marks_header.
        std::fprintf(f, "#SIZES cell=%zu vertex=%zu\n",
                     sizeof(Cell), sizeof(Vertex));
        std::fprintf(f, "#RADIUS %d\n", m_R);
      }
      std::fprintf(f, "OP %ld %s\n", id, kind);
      std::fprintf(f, "K %d\n", m_class);
      std::fprintf(f, "E %p %p\n", (void*)&*m_v0, (void*)&*m_v1);
      // What lock_zone() acquires today. Emitted as a set, not a count, so the
      // comparison is set-vs-set and a smaller measured set can be checked to
      // be a SUBSET rather than merely smaller.
      {
        std::unordered_set<const Vertex*> z;
        zone_today(z);
        for (const Vertex* w : z) std::fprintf(f, "T %p\n", (const void*)w);
      }
      // Each cell carries its four vertices, so the post-processor can turn
      // "this cell was touched" into "these vertices would have to be locked"
      // -- the spatial lock protects a cell by holding its vertices, so the
      // vertex set is the thing the zone is actually made of.
      for (const auto& s : m_cells)
      {
        int kl = m_ring.count(&*s.h) ? 0 : (m_mirror.count(&*s.h) ? 1 : 2);
        std::fprintf(f, "C %p %zu %d %d %p %p %p %p %d\n", (void*)&*s.h,
                     sizeof(Cell), s.depth, kl, (void*)&*s.v[0], (void*)&*s.v[1],
                     (void*)&*s.v[2], (void*)&*s.v[3], cell_coverage(s.h));
      }
      for (const auto& t : m_verts)
        std::fprintf(f, "V %p %zu %d %d\n", (void*)&*t.h, sizeof(Vertex),
                     t.depth, vertex_coverage(t.h));
      if (cfg.all_objects)
      {
        const int OUTSIDE = 90;
        for (auto c = m_tr.tds().cells().begin(); c != m_tr.tds().cells().end(); ++c)
          if (!m_pre_cells.count(&*c))
            std::fprintf(f, "C %p %zu %d 3 0 0 0 0 %d\n", (void*)&*c,
                         sizeof(Cell), OUTSIDE, cell_coverage(c));
        for (auto v = m_tr.tds().vertices().begin(); v != m_tr.tds().vertices().end(); ++v)
          if (!m_pre_verts.count(&*v))
            std::fprintf(f, "V %p %zu %d %d\n", (void*)&*v, sizeof(Vertex),
                         OUTSIDE, vertex_coverage(v));
      }
      std::fprintf(f, ".\n");
      std::fflush(f);
      mvlz_mark_begin = (unsigned long)id;     // <-- window opens here
      return;
    }

    self_test();
  }

  // A/A NEGATIVE CONTROL for diff mode: diffing the snapshot against itself,
  // with no operation in between, must report zero. A zero means nothing until
  // the check is shown able to reach non-zero AND to return zero when nothing
  // happened -- this campaign has already lost twelve "nulls" that were A vs A'.
  void self_test()
  {
    if (!m_active || !Mvlz_config::get().selftest) return;
    finish(m_kind + "_AA_null");
    m_active = true;                 // finish() only reads; the real end() runs
  }

  void end(bool /*op_returned_true*/)
  {
    if (!m_active) return;
    Mvlz_config& cfg = Mvlz_config::get();
    if (cfg.trace)
    {
      mvlz_mark_end = 1;             // <-- window closes here, before anything else
      m_active = false;
      if (cfg.exit_after >= 0 && ++windows_closed() >= cfg.exit_after)
      {
        if (cfg.manifest) { std::fprintf(cfg.manifest, "#DONE\n"); std::fflush(cfg.manifest); }
        std::_Exit(0);
      }
      return;
    }
    finish(m_kind);
  }

private:
  void finish(const std::string& kindname)
  {
    if (!m_active) return;
    m_active = false;
    Mvlz_totals::Kind& K = Mvlz_totals::get().kinds[kindname];
    ++K.ops;

    auto& cellc = m_tr.tds().cells();
    auto& vertc = m_tr.tds().vertices();
    int max_r = -1;
    bool wrote = false;
    std::unordered_set<const Vertex*> minzone;
    // every object whose bytes changed, or that the operation created --
    // gathered by the diff, so no write site has to be known in advance
    std::vector<const Cell*>   changed_cells, created_cells;
    std::vector<const Vertex*> changed_verts;

    for (const auto& s : m_cells)
    {
      if (!cellc.is_used(s.h))
      {
        ++K.cells_destroyed; wrote = true;
        max_r = (std::max)(max_r, s.depth);
        for (int i = 0; i < 4; ++i)
          if (!m_tr.is_infinite(s.v[i])) minzone.insert(&*s.v[i]);
        continue;
      }
      const unsigned char* now = reinterpret_cast<const unsigned char*>(&*s.h);
      if (std::memcmp(s.bytes.data(), now, sizeof(Cell)) == 0) continue;

      wrote = true; ++K.cells_written;
      changed_cells.push_back(&*s.h);
      max_r = (std::max)(max_r, s.depth);
      for (std::size_t i = 0; i < sizeof(Cell); ++i)
        if (s.bytes[i] != now[i]) ++K.cell_offsets[(int)i];
      for (int i = 0; i < 4; ++i)
      {
        if (s.v[i] != s.h->vertex(i))   ++K.f_vertex;
        if (s.n[i] != s.h->neighbor(i)) ++K.f_neighbor;
        if (!m_tr.is_infinite(s.h->vertex(i))) minzone.insert(&*s.h->vertex(i));
        if (!m_tr.is_infinite(s.v[i]))         minzone.insert(&*s.v[i]);
      }
      if      (m_ring.count(&*s.h))   ++K.wc_ring;
      else if (m_mirror.count(&*s.h)) ++K.wc_mirror;
      else if (s.depth <= 1)          ++K.wc_star_other;
      else                            ++K.wc_far;
    }

    for (const auto& t : m_verts)
    {
      if (!vertc.is_used(t.h))
      {
        ++K.verts_destroyed; wrote = true;
        if (!t.infinite) { max_r = (std::max)(max_r, t.depth); minzone.insert(&*t.h); }
        continue;
      }
      const unsigned char* now = reinterpret_cast<const unsigned char*>(&*t.h);
      if (std::memcmp(t.bytes.data(), now, sizeof(Vertex)) == 0) continue;
      wrote = true;
      if (t.infinite) { ++K.inf_vertex_written; continue; }
      ++K.verts_written;
      changed_verts.push_back(&*t.h);
      max_r = (std::max)(max_r, t.depth);
      for (std::size_t i = 0; i < sizeof(Vertex); ++i)
        if (t.bytes[i] != now[i]) ++K.vert_offsets[(int)i];
      if (t.cell != t.h->cell()) ++K.f_cell;
      if (t.x != t.h->point().x() || t.y != t.h->point().y()
          || t.z != t.h->point().z()) ++K.f_point;
      minzone.insert(&*t.h);
    }

    // objects that did not exist before
    {
      std::vector<Cell_handle> inc;
      std::unordered_set<const Cell*> post;
      std::unordered_set<const Vertex*> fresh;
      for (Vertex_handle v : {m_v0, m_v1})
      {
        if (!vertc.is_used(v)) continue;
        inc.clear();
        m_tr.incident_cells(v, std::back_inserter(inc));
        for (Cell_handle c : inc)
          if (post.insert(&*c).second && !m_pre_cells.count(&*c))
          {
            ++K.cells_created;
            changed_cells.push_back(&*c);
            created_cells.push_back(&*c);
            for (int i = 0; i < 4; ++i)
            {
              Vertex_handle w = c->vertex(i);
              if (m_tr.is_infinite(w)) continue;
              minzone.insert(&*w);
              if (m_pre_verts.count(&*w))
                max_r = (std::max)(max_r, depth_of(w));
              else
                // A vertex the operation itself created. It has no BFS depth
                // because it did not exist when the ball was built, and
                // depth_of() would hand back the not-found sentinel R+1 and
                // report a spurious write at the ball bound -- which is
                // exactly what the first run did. It sits on the element by
                // construction, so its depth is 0 and it is counted once.
                if (fresh.insert(&*w).second) ++K.verts_created;
            }
          }
      }
    }

    // ---- LOCK COVERAGE ----------------------------------------------------
    // Runs here, before the executor's unlock_all_elements(): the zone this
    // operation acquired is still held. Every object the diff says changed
    // must be covered; the diff is complete by construction, so this is a
    // sufficiency test of the zone rather than of a list of write sites.
    if (Mvlz_config::get().lock_check)
    {
      if (m_tr.get_lock_data_structure() == nullptr) ++K.lc_no_lock_ds;
      else
      {
        bool viol = false;
        std::unordered_set<const Cell*> is_new(created_cells.begin(), created_cells.end());
        auto* lds = m_tr.get_lock_data_structure();
        for (const Cell*  cp : changed_cells)
        {
          ++K.lc_cells_checked;
          Cell_handle ch(const_cast<Cell*>(cp));
          if (cell_covered(ch, K.lc_infinite_skipped)) continue;
          ++K.lc_cells_uncovered; viol = true;
          if (is_new.count(cp)) ++K.lc_uncov_created; else ++K.lc_uncov_modified;
          // which vertex was not held -- one the operation created, or one
          // that already existed? A created vertex sits at a NEW position, and
          // the lock grid is keyed on position.
          for (int i = 0; i < 4; ++i)
          {
            Vertex_handle v = ch->vertex(i);
            if (v == Vertex_handle() || m_tr.is_infinite(v)) continue;
            if (lds->is_locked_by_this_thread(v->point())) continue;
            if (m_pre_verts.count(&*v)) ++K.lc_uncov_by_old_vertex;
            else                        ++K.lc_uncov_by_new_vertex;
          }
        }
        for (const Vertex* vp : changed_verts)
        {
          ++K.lc_verts_checked;
          Vertex_handle vh(const_cast<Vertex*>(vp));
          if (!vertex_covered(vh, K.lc_infinite_skipped))
          { ++K.lc_verts_uncovered; viol = true; }
        }
        if (viol) ++K.lc_ops_with_violation;
      }
    }

    if (wrote) ++K.ops_written;
    if (max_r >= 0)
    {
      if (max_r >= (int)K.radius.size()) max_r = (int)K.radius.size() - 1;
      ++K.radius[max_r];
      if (max_r >= m_R) ++K.saturated;
    }
    K.zone_now_verts += m_zone_now;
    K.zone_min_verts += (long)minzone.size();
  }
};

struct Mvlz_reporter
{
  ~Mvlz_reporter()
  {
    if (Mvlz_config::get().manifest) std::fclose(Mvlz_config::get().manifest);
    const char* p = std::getenv("CGAL_TR_MVLZ_OUT");
    if (!Mvlz_config::get().trace)
      Mvlz_totals::get().report(p ? p : "mvlz.txt");
  }
};

inline Mvlz_reporter& mvlz_reporter()
{
  (void)mvlz_marks_header;   // keep the static initialiser
  // Touch the totals and the config FIRST so they are constructed before the
  // reporter and therefore destroyed after it. Without this the reporter's
  // destructor runs against an already-destroyed std::map: the first run
  // printed its results and then spun forever on a zero-byte report file.
  Mvlz_config::get();
  Mvlz_totals::get();
  static Mvlz_reporter r;
  return r;
}

} // internal
} // Tetrahedral_remeshing
} // CGAL

#endif // CGAL_TR_MVLZ_PROBE
#endif // CGAL_TETRAHEDRAL_REMESHING_MVLZ_PROBE_H
