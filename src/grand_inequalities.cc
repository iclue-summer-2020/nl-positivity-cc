// Copyright 2020 ICLUE @ UIUC. All rights reserved.

#include <algorithm>
#include <climits>
#include <deque>
#include <vector>
#include <stack>

#include <combinations.hpp>
#include <nl_positivity/grand_inequalities.h>
#include <nlnum/nlnum.h>
#include <product.hpp>
#include <range.hpp>

namespace nl_positivity {

std::vector<std::vector<Int>> product(const Set& s, Int repeat) {
  std::vector<std::vector<Int>> ans;
  std::stack<std::vector<Int>> call_stack;
  call_stack.push({});

  while (!call_stack.empty()) {
    const auto top = call_stack.top();
    call_stack.pop();

    if (top.size() >= repeat) {
      ans.push_back(top);
      continue;
    }

    for (const auto& e : s) {
      std::vector<Int> vt = top;
      vt.push_back(e);
      call_stack.push(vt);
    }
  }

  return ans;
}

std::vector<std::vector<Set>> disjoints(const Int n, const Int nsets) {
  std::vector<std::vector<Set>> ds{};
  const auto N = iter::range(nsets+1);
  for (const auto& bits : product(Set{N.begin(), N.end()}, n)) {
    std::deque<Set> S(nsets+1);
    for (size_t i = 0; i < bits.size(); ++i) {
      S[bits[i]].insert(i + 1);
    }
    S.pop_front();
    ds.emplace_back(S.begin(), S.end());
  }

  return ds;
}

std::vector<Int> tau(const Set& X) {
  const size_t d = X.size();

  std::vector<Int> ss{X.begin(), X.end()};
  std::sort(ss.begin(), ss.end());

  const auto R = iter::range(1, static_cast<int32_t>(d)+1);
  const std::vector<Int> rr{R.begin(), R.end()};

  std::vector<Int> tt;
  for (size_t idx = 0; idx < d; ++idx) {
    tt.push_back(ss[idx]-rr[idx]);
  }
  std::reverse(tt.begin(), tt.end());

  return tt;
}

template<typename T>
std::vector<Set> collect(const T& t) {
  std::vector<Set> ans{};
  for (const auto& e : t) {
    ans.push_back({e.begin(), e.end()});
  }
  return ans;
}

// Returns the empty set in the case that k == 0.
std::vector<Set> combinations(const std::vector<Int> R, const size_t k) {
  if (k == 0) return {{}};
  return collect(iter::combinations(R, k));
}

bool is_good(Int n, const std::vector<Int>& R, const Set& S, const Set& A,
             const Set& B, const Set& C, const Set& Ap, const Set& Bp,
             const Set& Cp, Sets* s) {
  const auto AV = combinations(R, Ap.size());
  for (const auto& a : iter::product<2>(AV)) {
    const Set& A1 = std::get<0>(a);
    const Set& A2 = std::get<1>(a);
    if (nlnum::lrcoef(tau(Ap), tau(A1), tau(A2)) <= 0) continue;

    const auto BV = combinations(R, Bp.size());
    for (const auto& b : iter::product<2>(BV)) {
      const Set& B1 = std::get<0>(b);
      const Set& B2 = std::get<1>(b);
      if (nlnum::lrcoef(tau(Bp), tau(B1), tau(B2)) <= 0) continue;

      const auto CV = combinations(R, Cp.size());
      for (const auto& c : iter::product<2>(CV)) {
        const Set& C1 = std::get<0>(c);
        const Set& C2 = std::get<1>(c);
        if (nlnum::lrcoef(tau(Cp), tau(C1), tau(C2)) <= 0) continue;

        const int32_t mA = static_cast<int32_t>(std::min(Bp.size(), Cp.size()));
        const int32_t mB = static_cast<int32_t>(std::min(Ap.size(), Cp.size()));
        const int32_t mC = static_cast<int32_t>(std::min(Ap.size(), Bp.size()));

        Set Ac{};
        Set Bc{};
        Set Cc{};
        Set A1c{};
        Set B1c{};
        Set C1c{};
        Set A2c{};
        Set B2c{};
        Set C2c{};
        std::vector<std::pair<const Set&, Set&>> SS = {
            {A, Ac},   {B, Bc},   {C, Cc},   {A1, A1c}, {B1, B1c},
            {C1, C1c}, {A2, A2c}, {B2, B2c}, {C2, C2c},
        };
        for (auto& it : SS) {
          std::set_difference(S.begin(), S.end(), it.first.begin(),
                              it.first.end(),
                              std::inserter(it.second, it.second.begin()));
        }

        Set Aw{};
        Set Bw{};
        Set Cw{};
        Set A1w{};
        Set B1w{};
        Set C1w{};
        Set A2w{};
        Set B2w{};
        Set C2w{};
        std::vector<std::tuple<const Set&, const Set&, int32_t, Set&>> TT = {
            {A, Ac, mA, Aw},    {B, Bc, mB, Bw},    {C, Cc, mC, Cw},
            {A1, A1c, mC, A1w}, {B1, B1c, mA, B1w}, {C1, C1c, mB, C1w},
            {A2, A2c, mB, A2w}, {B2, B2c, mC, B2w}, {C2, C2c, mA, C2w},
        };
        for (auto& it : TT) {
          const Set& X = std::get<0>(it);
          const Set& Xc = std::get<1>(it);
          const int32_t m = std::get<2>(it);
          Set& Xw = std::get<3>(it);
          const auto nn = static_cast<int32_t>(n);
          const auto Rq =
              iter::range(nn + 1, nn + static_cast<int32_t>(X.size()) - m + 1);
          const Set Rs{Rq.begin(), Rq.end()};
          std::set_union(Xc.begin(), Xc.end(), Rs.begin(), Rs.end(),
                         std::inserter(Xw, Xw.begin()));
        }

        if (nlnum::lrcoef(tau(Cw), tau(A1w), tau(B2w)) <= 0 ||
            nlnum::lrcoef(tau(Bw), tau(C1w), tau(A2w)) <= 0 ||
            nlnum::lrcoef(tau(Aw), tau(B1w), tau(C2w)) <= 0)
          continue;

        if (s != nullptr) {
          s->A = A;
          s->B = B;
          s->C = C;
          s->Ap = Ap;
          s->Bp = Bp;
          s->Cp = Cp;
          s->A1 = A1;
          s->B1 = B1;
          s->C1 = C1;
          s->A2 = A2;
          s->B2 = B2;
          s->C2 = C2;
        }
        return true;
      }
    }
  }
  return false;
}

std::vector<Sets> grand_ineqs(const Int n) {
  const auto djs = disjoints(n, 2);
  const auto RR = iter::range(1, static_cast<int32_t>(n)+1);
  const std::vector<int32_t> Rv{RR.begin(), RR.end()};
  const std::vector<Int> R{Rv.begin(), Rv.end()};
  const Set S{R.begin(), R.end()};

  std::vector<Sets> ans{};
#pragma omp parallel for schedule(dynamic)
  for (auto ita = djs.begin(); ita < djs.end(); ++ita) {
    const auto& Aq = *ita;
    const Set& A = Aq[0];
    const Set& Ap = Aq[1];
    for (const auto& Bq : djs) {
      const Set& B = Bq[0];
      const Set& Bp = Bq[1];
      for (const auto& Cq : djs) {
        const Set& C = Cq[0];
        const Set& Cp = Cq[1];

        if (A.size() < std::max(Bp.size(), Cp.size()) ||
            B.size() < std::max(Ap.size(), Cp.size()) ||
            C.size() < std::max(Ap.size(), Bp.size()))
          continue;

        Sets s;
        if (is_good(n, R, S, A, B, C, Ap, Bp, Cp, &s)) {
#pragma omp critical
          ans.push_back(s);
        }
      }
    }
  }

  return ans;
}

}  // namespace nl_positivity
