// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#include <gflags/gflags.h>
#include <nl_positivity/grand_inequalities.h>
#include <prettyprint.hpp>

#include <iostream>
#include <map>
#include <tuple>
#include <vector>

DEFINE_uint32(n, 0, "n");
DEFINE_uint32(k, 0, "k");

using nl_positivity::flagger;
using nl_positivity::grand_ineqs;
using nl_positivity::Set;
using nl_positivity::Sets;

int main(int argc, char** argv) {
  gflags::SetUsageMessage("Finds counterexamples to the NL numbers claim");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  const auto n = FLAGS_n;
  const auto k = FLAGS_k;

  const auto setss = grand_ineqs(n, [](int64_t c) -> bool { return c == 1; });
  std::cout << "Number of sets: " << setss.size() << std::endl;
  std::cout << "Sets:" << std::endl;

  for (const auto& sets : setss) {
    std::cout << sets << std::endl;
  }

  std::cout << "===========================" << std::endl;

  const auto flagged = flagger(n, k);
  std::cout << "Number flagged: " << flagged.size() << std::endl;

  for (const auto& flag : flagged) {
    std::cout << flag << std::endl;
  }

  return EXIT_SUCCESS;
}
