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

#ifndef CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H
#define CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>

#include <algorithm>
#include <functional>
#include <iterator>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#include <CGAL/Real_timer.h>
#include <cstddef>
#include <iostream>
#endif

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

/**
* Base class of the elementary remeshing operations.
*
* An operation is described by
*  - `elements()`          : the range scanned to collect candidates,
*  - `predicate()`         : whether the operation applies to an element, and
*                            with which priority (`std::nullopt` = do not apply),
*  - `execute_operation()` : the modification itself,
*  - `affected_elements()` : for operations that need re-queueing, the local
*                            range of elements whose priority may have changed
*                            after a successful operation.
*
* The executor drives the loop : it scans `elements()`, keeps what `predicate()`
* accepts, and - for re-queueing operations - runs the very same `predicate()`
* on `affected_elements()` after each successful operation. Using one predicate
* for both ranges is what keeps re-queued candidates consistent with the
* initially collected ones.
*
* \tparam ElementType durable representation of an element, i.e. what the
*   executor stores in its work list. It may differ from the type yielded by
*   `elements()`, in which case `to_element()` must be overridden.
*/
template <typename C3t3_,
          typename ElementType,
          typename ElementRange,
          typename Priority = typename C3t3_::Triangulation::Geom_traits::FT,
          typename Result = bool,
          typename ElementCompare = std::less<ElementType> >
class Elementary_operation
{
public:
  using C3t3 = C3t3_;
  using Triangulation = typename C3t3::Triangulation;
  using Element_type = ElementType;
  using Element_range = ElementRange;
  using Priority_type = Priority;
  using Result_type = Result;
  using Element_compare = ElementCompare;

  // what `elements()` and `affected_elements()` yield
  using Scan_type = typename std::decay<
      decltype(*std::begin(std::declval<const ElementRange&>()))>::type;
  using Affected_range = std::vector<Scan_type>;

  /**
  * Sink used to drop the elements that an operation destroys while it runs.
  * Handles of destroyed cells get recycled, so such elements must leave the
  * work list *during* the operation : afterwards they can no longer be
  * enumerated, and probing them would read recycled memory.
  */
  class Invalidation_sink
  {
  public:
    virtual ~Invalidation_sink() = default;
    virtual void invalidate(const Scan_type& e) = 0;
  };

  Elementary_operation() = default;
  virtual ~Elementary_operation() = default;

  /// range scanned to collect the initial candidates
  virtual Element_range elements(const C3t3& c3t3) const = 0;

  /// does the operation apply to `e`, and with which priority?
  virtual std::optional<Priority_type> predicate(const Scan_type& e,
                                                 const C3t3& c3t3) const = 0;

  /// converts a scanned element into the representation stored by the executor
  virtual Element_type to_element(const Scan_type& e) const
  {
    if constexpr (std::is_convertible_v<Scan_type, Element_type>)
      return static_cast<Element_type>(e);
    else
    {
      CGAL_assertion(false); // must be overridden when the two types differ
      return Element_type();
    }
  }

  virtual Result_type execute_operation(const Element_type& e, C3t3& c3t3) = 0;

  /// same, for the operations that destroy elements of the work list
  virtual Result_type execute_operation(const Element_type& e, C3t3& c3t3,
                                        Invalidation_sink&)
  {
    return execute_operation(e, c3t3);
  }

  virtual bool succeeded(const Result_type& r) const { return r != Result_type(); }

  /// elements whose priority may have changed after a successful operation
  virtual Affected_range affected_elements(const Result_type&, const C3t3&) const
  {
    return Affected_range();
  }

  virtual bool requires_ordered_processing() const { return false; }
  virtual bool requires_requeue() const { return false; }
  /// for ordered operations : is the smallest priority processed first?
  virtual bool process_smallest_first() const { return true; }

  virtual std::string operation_name() const = 0;
};


template <typename Operation>
class Elementary_operation_execution_sequential
{
public:
  using C3t3 = typename Operation::C3t3;
  using Element_type = typename Operation::Element_type;
  using Priority_type = typename Operation::Priority_type;
  using Element_compare = typename Operation::Element_compare;
  using Scan_type = typename Operation::Scan_type;

  bool execute(Operation& op, C3t3& c3t3) const
  {
    return op.requires_requeue() ? execute_with_requeue(op, c3t3)
                                 : execute_static(op, c3t3);
  }

private:
  /// keeps the elements of `op.elements()` that `op.predicate()` accepts
  template <typename Insert>
  void collect(const Operation& op, const C3t3& c3t3, Insert insert) const
  {
    for (const auto& s : op.elements(c3t3))
    {
      const std::optional<Priority_type> p = op.predicate(s, c3t3);
      if (p != std::nullopt)
        insert(op.to_element(s), p.value());
    }
  }

  /// one pass over a candidate list collected once
  bool execute_static(Operation& op, C3t3& c3t3) const
  {
    using Candidate = std::pair<Element_type, Priority_type>;
    std::vector<Candidate> candidates;
    collect(op, c3t3, [&candidates](const Element_type& e, const Priority_type& p)
                      { candidates.push_back(Candidate(e, p)); });

    if (candidates.empty())
      return false;

    if (op.requires_ordered_processing())
    {
      const bool smallest_first = op.process_smallest_first();
      std::stable_sort(candidates.begin(), candidates.end(),
        [smallest_first](const Candidate& a, const Candidate& b)
        { return smallest_first ? (a.second < b.second) : (a.second > b.second); });
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::size_t nb_done = 0;
    const std::size_t nb_candidates = candidates.size();
    CGAL::Real_timer timer;
    timer.start();
#endif
    for (const Candidate& c : candidates)
    {
      if (op.succeeded(op.execute_operation(c.first, c3t3)))
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        ++nb_done;
#endif
      }
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    timer.stop();
    std::cout << op.operation_name() << ": " << nb_done << "/"
              << nb_candidates << " done ("
              << timer.time() << " sec)." << std::endl;
#endif
    return true;
  }

  /**
  * Priority driven pass : always process the best element of the work list,
  * then re-evaluate the elements it affected, so that those that became
  * candidates are processed within the same phase.
  */
  bool execute_with_requeue(Operation& op, C3t3& c3t3) const
  {
    // work list keyed by priority; a "largest first" operation stores the
    // opposite priority, so that the best element is always the smallest key
    using Bimap = boost::bimap<
        boost::bimaps::set_of<Element_type, Element_compare>,
        boost::bimaps::multiset_of<Priority_type, std::less<Priority_type> > >;

    const bool smallest_first = op.process_smallest_first();
    const auto key = [smallest_first](const Priority_type& p)
                     { return smallest_first ? p : -p; };

    Bimap worklist;
    const auto update = [&worklist](const Element_type& e, const Priority_type& k)
    {
      auto it = worklist.left.find(e);
      if (it != worklist.left.end()) worklist.left.replace_data(it, k);
      else worklist.left.insert(typename Bimap::left_map::value_type(e, k));
    };
    const auto erase = [&worklist](const Element_type& e)
    {
      auto it = worklist.left.find(e);
      if (it != worklist.left.end()) worklist.left.erase(it);
    };

    collect(op, c3t3, [&update, &key](const Element_type& e, const Priority_type& p)
                      { update(e, key(p)); });
    if (worklist.empty())
      return false;

    // elements destroyed by an operation leave the work list at once
    class Sink : public Operation::Invalidation_sink
    {
    public:
      Sink(Bimap& w, const Operation& o) : m_worklist(w), m_op(o) {}
      void invalidate(const Scan_type& s) override
      {
        auto it = m_worklist.left.find(m_op.to_element(s));
        if (it != m_worklist.left.end())
          m_worklist.left.erase(it);
      }
    private:
      Bimap& m_worklist;
      const Operation& m_op;
    };
    Sink sink(worklist, op);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::size_t nb_done = 0;
    CGAL::Real_timer timer;
    timer.start();
#endif
    while (!worklist.empty())
    {
      const Element_type e = worklist.right.begin()->second;
      worklist.right.erase(worklist.right.begin());

      const typename Operation::Result_type result
        = op.execute_operation(e, c3t3, sink);
      if (!op.succeeded(result))
        continue;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      ++nb_done;
#endif
      // re-evaluate what the operation affected : an element that became a
      // candidate enters the work list, one that stopped being a candidate
      // (its priority changed with the mesh) leaves it
      for (const Scan_type& a : op.affected_elements(result, c3t3))
      {
        const std::optional<Priority_type> p = op.predicate(a, c3t3);
        if (p != std::nullopt) update(op.to_element(a), key(p.value()));
        else                   erase (op.to_element(a));
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

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H
