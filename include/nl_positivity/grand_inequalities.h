// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#ifndef NL_POSITIVITY_GRAND_INEQUALITIES_H_
#define NL_POSITIVITY_GRAND_INEQUALITIES_H_

#include <cstdint>
#include <set>
#include <vector>

namespace nl_positivity {

typedef uint64_t Int;
typedef std::set<Int> Set;

struct Sets {
  Set A, B, C, Ap, Bp, Cp, A1, B1, C1, A2, B2, C2;
  Sets() {}
  Sets(
      const Set& A,
      const Set& B,
      const Set& C,
      const Set& Ap,
      const Set& Bp,
      const Set& Cp,
      const Set& A1,
      const Set& B1,
      const Set& C1,
      const Set& A2,
      const Set& B2,
      const Set& C2
  ) : A{A}, B{B}, C{C}, Ap{Ap}, Bp{Bp}, Cp{Cp},
    A1{A1}, B1{B1}, C1{C1}, A2{A2}, B2{B2}, C2{C2} { }
};

// Computes the generalized Cartesian product of a set with itself.
std::vector<std::vector<Int>> product(const Set& s, Int repeat);

std::vector<Sets> grand_ineqs(Int n, bool (*cond)(int64_t));

// Computes the tau function defined in the NL paper.
std::vector<Int> tau(const Set& X);

// Returns all `nsets`-tuple of sets of [n] that are disjoint.
std::vector<std::vector<Set>> disjoints(Int n, Int nsets);

}  // namespace nl_positivity

#endif  // NL_POSITIVITY_GRAND_INEQUALITIES_H_
