// Copyright 2020 ICLUE @ UIUC. All rights reserved.

#include <algorithm>
#include <climits>
#include <deque>
#include <vector>
#include <stack>
#include <unordered_set>

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
  std::vector<std::vector<Set>> ds;
  const auto N = iter::range(nsets+1);
  for (const auto& bits : product(Set{N.begin(), N.end()}, n)) {
    std::deque<Set> S(nsets+1);
    for (size_t i = 0; i < bits.size(); ++i) {
      S[bits[i]].insert(i + 1);
    }
    S.pop_front();
    ds.push_back({S.begin(), S.end()});
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
  std::vector<Set> ans;
  for (const auto& e : t) {
    ans.push_back({e.begin(), e.end()});
  }
  return ans;
}

std::vector<Sets> grand_ineqs(const Int n) {
  const auto djs = disjoints(n, 2);
  const auto RR = iter::range(1, static_cast<int32_t>(n)+1);
  const std::vector<int32_t> Rv{RR.begin(), RR.end()};
  const std::vector<Int> R{Rv.begin(), Rv.end()};

  for (const auto& Aq : djs) {
    const Set& A = Aq[0];
    const Set& Ap = Aq[1];
    for (const auto& Bq : djs) {
      const Set& B = Bq[0];
      const Set& Bp = Bq[1];
      for (const auto& Cq : djs) {
        const Set& C = Cq[0];
        const Set& Cp = Cq[1];

        if (A.size() < std::max(Bp.size(), Cp.size())
          || B.size() < std::max(Ap.size(), Bp.size())
          || C.size() < std::max(Ap.size(), Cp.size())) continue;

        const auto AV = collect(iter::combinations(R, Ap.size()));
        for (const auto& a : iter::product<2>(AV)) {
          const Set& A1 = std::get<0>(a);
          const Set& A2 = std::get<1>(a);
          if (nlnum::lrcoef(tau(Ap), tau(A1), tau(A2)) <= 0) continue;

          const auto BV = collect(iter::combinations(R, Bp.size()));
          for (const auto& b : iter::product<2>(BV)) {
            const Set& B1 = std::get<0>(b);
            const Set& B2 = std::get<1>(b);
            if (nlnum::lrcoef(tau(Bp), tau(B1), tau(B2)) <= 0) continue;

            const auto CV = collect(iter::combinations(R, Cp.size()));
            for (const auto& c : iter::product<2>(CV)) {
              const Set& C1 = std::get<0>(c);
              const Set& C2 = std::get<1>(c);
              if (nlnum::lrcoef(tau(Cp), tau(C1), tau(C2)) <= 0) continue;
            }
          }
        }
      }
    }
  }
  return {};
}

}  // namespace nl_positivity
