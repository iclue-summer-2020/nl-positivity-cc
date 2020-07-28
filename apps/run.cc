// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#include <gflags/gflags.h>
#include <nl_positivity/grand_inequalities.h>
#include <prettyprint.hpp>

#include <iostream>
#include <map>
#include <tuple>
#include <vector>

DEFINE_uint32(n, 0, "n");

using nl_positivity::flagger;
using nl_positivity::grand_ineqs;
using nl_positivity::Set;
using nl_positivity::Sets;

int main(int argc, char** argv) {
  gflags::SetUsageMessage("Finds counterexamples to the NL numbers claim");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  const auto n = FLAGS_n;

  const auto setsA = grand_ineqs(n, [](int64_t c) -> bool { return c > 0; });
  const auto setsB = grand_ineqs(n, [](int64_t c) -> bool { return c == 1; });

  std::cout << setsA.size() << ", " << setsB.size() << std::endl;

  // We want to see which {A,B,C} sets only have one entry.
  std::map<std::tuple<Set, Set, Set, Set, Set, Set>, std::vector<Sets>> table;
  for (const auto& s : setsA) {
    auto& v = table[{s.A, s.B, s.C, s.Ap, s.Bp, s.Cp}];
    v.push_back(s);
  }
  for (const auto& s : setsB) {
    auto& v = table[{s.A, s.B, s.C, s.Ap, s.Bp, s.Cp}];
    v.push_back(s);
  }
  for (const auto& kv : table) {
    if (kv.second.size() == 1) {
      std::cout << kv.second << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
