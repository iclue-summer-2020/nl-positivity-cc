// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#include <gflags/gflags.h>
#include <nl_positivity/grand_inequalities.h>
#include <prettyprint.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>

DEFINE_uint32(n, 0, "n");
DEFINE_uint32(k, 0, "k");

int main(int argc, char** argv) {
  gflags::SetUsageMessage("Finds counterexamples to the NL numbers claim");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  const auto n = FLAGS_n;
  const auto k = FLAGS_k;
  const auto flags = nl_positivity::flagger(n, k);

  std::cout << flags.size() << std::endl;
  if (!flags.empty()) {
    std::cout << flags << std::endl;
  }
  return EXIT_SUCCESS;
}
