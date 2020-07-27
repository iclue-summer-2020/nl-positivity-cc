// Copyright (c) 2020 [Your Name]. All rights reserved.

#define CATCH_CONFIG_MAIN

#include <exception>
#include <iostream>
#include <numeric>

#include <catch2/catch.hpp>
#include <nl_positivity/grand_inequalities.h>

TEST_CASE("n=1", "") {
}

TEST_CASE("product", "") {
  const nl_positivity::Set s = {1, 2, 3};
  const auto ps = nl_positivity::product(s, 2);
  REQUIRE(ps.size() == 9);
}

TEST_CASE("disjoints", "") {
  const auto djs = nl_positivity::disjoints(3, 2);
  REQUIRE(djs.size() == 27);
}

TEST_CASE("tau", "") {
  const nl_positivity::Set s = {4, 6, 3, 9, 1};
  const auto tt = nl_positivity::tau(s);
  const std::vector<nl_positivity::Int> expected = {4, 2, 1, 1, 0};
  REQUIRE(tt == expected);
}

TEST_CASE("grand", "") {
  const auto x = nl_positivity::grand_ineqs(3);
  REQUIRE(x.size() == 0);
}
