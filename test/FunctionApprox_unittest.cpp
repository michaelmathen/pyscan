//
// Created by mmath on 10/2/17.
//


#include "../src/RectangleScan.hpp"
#include "../src/FunctionApprox.hpp"
#include "../src/DiskScan.hpp"
#include "Utilities.hpp"

#include <limits.h>
#include <random>
#include <iostream>

#include "gtest/gtest.h"
namespace {

// Step 2. Use the TEST macro to define your tests.
//
// TEST has two parameters: the test case name and the test name.
// After using the macro, you should define your test logic between a
// pair of braces.  You can use a bunch of macros to indicate the
// success or failure of a test.  EXPECT_TRUE and EXPECT_EQ are
// examples of such macros.  For a complete list, see gtest.h.
//
// <TechnicalDetails>
//
// In Google Test, tests are grouped into test cases.  This is how we
// keep test code organized.  You should put logically related tests
// into the same test case.
//
// The test case name and the test name should both be valid C++
// identifiers.  And you should not use underscore (_) in the names.
//
// Google Test guarantees that each test you define is run exactly
// once, but it makes no guarantee on the order the tests are
// executed.  Therefore, you should write your tests in such a way
// that their results don't depend on their order.
//
// </TechnicalDetails>

TEST(ApproximateHullTest, Kulldorff) {
    const static int test_size = 10000;
    double rho = .001;
    double alpha = .001;
    double eps = .01;
    auto pts = pyscantest::randomVec(test_size);

    auto avg = [&] (pyscan::VecD const& v1, pyscan::VecD const& v2) {
        pyscan::VecD v_out = v1 + v2;
        return v_out * 1.0 / mag(v_out);
    };

    double maxV_approx = pyscan::approximateHull(eps,
      [&](pyscan::VecD pt) {
        return pyscan::regularized_kulldorff(pt[0], pt[1], rho);
      },
      [&] (pyscan::VecD dir){
        return pyscantest::maxVecD(pts, [&](pyscan::VecD const& v) {
          return pyscan::dot(v, dir);
      });
    });

    auto maxV_exact_pt = pyscantest::maxVecD(pts, [&](pyscan::VecD const& pt){
      return pyscan::regularized_kulldorff(pt[0], pt[1], rho);
    });

    double maxV_exact = pyscan::regularized_kulldorff(maxV_exact_pt[0], maxV_exact_pt[1], rho);
    std::cout << maxV_approx << " " << maxV_exact << std::endl;
    EXPECT_NEAR(maxV_exact, maxV_approx, eps);

  }

}
// Step 3. Call RUN_ALL_TESTS() in main().
//
// We do this by linking in src/gtest_main.cc file, which consists of
// a main() function which calls RUN_ALL_TESTS() for us.
//
// This runs all the tests you've defined, prints the result, and
// returns 0 if successful, or 1 otherwise.
//
// Did you notice that we didn't register the tests?  The
// RUN_ALL_TESTS() macro magically knows about all the tests we
// defined.  Isn't this convenient?
