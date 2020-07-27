// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#include <gflags/gflags.h>
#include <nl_positivity/grand_inequalities.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>

DEFINE_uint32(n, 0, "n");

int main(int argc, char** argv) {
  gflags::SetUsageMessage("Prints each partition of n");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_n <= 0) {
    std::cerr << "Please provide a value of n" << std::endl;
    return EXIT_FAILURE;
  }

  if (FLAGS_n >= (1u << 8u)) {
    std::cerr << "This value of n is too large" << std::endl;
    return EXIT_FAILURE;
  }

  const auto n = static_cast<int8_t>(FLAGS_n);

  std::vector<nl_positivity::Sets> sets = nl_positivity::grand_ineqs(n);
  std::cout << sets.size() << std::endl;
  return EXIT_SUCCESS;
}
