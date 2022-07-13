#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

TEST_CASE("Always Pass") { REQUIRE(1 == 1); }