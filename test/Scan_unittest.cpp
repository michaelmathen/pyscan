//
// Created by mmath on 10/2/17.
//


#include "../src/RectangleScan.hpp"
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


    TEST(DiskTest, Kulldorff) {

        const static int n_size = 50;
        const static int s_size = 1000;
        const static double rho = .01;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomPoints(s_size);
        auto b_pts = pyscantest::randomPoints(s_size);
        auto d1 = pyscan::diskScanSlowStat(n_pts, m_pts, b_pts, rho);
        auto d2 = pyscan::diskScanStat(n_pts, m_pts, b_pts, rho);
        std::cout << d1.getA() << " " << d1.getB() << " " << d1.getR() << " " << d1.fValue() << std::endl;
        std::cout << d2.getA() << " " << d2.getB() << " " << d2.getR() << " " << d2.fValue() << std::endl;

        auto f = [&](double m, double b) {
            return pyscan::kulldorff(m, b, rho);
        };
        EXPECT_FLOAT_EQ(d1.fValue(), evaluateRegion(m_pts, b_pts, d1, f));
        EXPECT_FLOAT_EQ(d2.fValue(), evaluateRegion(m_pts, b_pts, d2, f));
        EXPECT_FLOAT_EQ(d1.getA(), d2.getA());
        EXPECT_FLOAT_EQ(d1.getB(), d2.getB());
        EXPECT_FLOAT_EQ(d1.getR(), d2.getR());

    }

    TEST(DiskTestLabel, Kulldorff) {

        const static int n_size = 50;
        const static int s_size = 1000;
        const static double rho = .01;
        auto n_pts = pyscantest::randomLPoints(n_size, 10);
        auto m_pts = pyscantest::randomLPoints(s_size, 10);
        auto b_pts = pyscantest::randomLPoints(s_size, 10);
        auto d1 = pyscan::diskScanSlowStatLabels(n_pts, m_pts, b_pts, rho);
        auto d2 = pyscan::diskScanStatLabels(n_pts, m_pts, b_pts, rho);
        std::cout << d1.getA() << " " << d1.getB() << " " << d1.getR() << " " << d1.fValue() << std::endl;
        std::cout << d2.getA() << " " << d2.getB() << " " << d2.getR() << " " << d2.fValue() << std::endl;

        auto f = [&](double m, double b) {
            return pyscan::kulldorff(m, b, rho);
        };
        EXPECT_FLOAT_EQ(d1.fValue(), evaluateRegion(m_pts, b_pts, d1, f));
        EXPECT_FLOAT_EQ(d2.fValue(), evaluateRegion(m_pts, b_pts, d2, f));
        EXPECT_FLOAT_EQ(d1.getA(), d2.getA());
        EXPECT_FLOAT_EQ(d1.getB(), d2.getB());
        EXPECT_FLOAT_EQ(d1.getR(), d2.getR());

    }
    //TODO write labeled disk scan tests.


    TEST(ApproximateHullTest, Kulldorff) {
        const static int test_size = 1000;
        auto pts = pyscantest::randomVec(test_size);

        auto avg = [&] (pyscan::VecD const& v1, pyscan::VecD const& v2) {
            pyscan::VecD v_out = v1 + v2;
            return v_out * 1.0 / mag(v_out);
        };

        auto line_f = [&] (pyscan::VecD dir){
            pyscan::VecD max_pt(0.0, 0.0);
            for (auto pt : pts) {
                if (dot(pt, dir) > dot(max_pt, dir)) {
                    max_pt = pt;
                }
            }
            return max_pt;
        };
        /*
        EXPECT_EQ(1, Factorial(-5));
        EXPECT_EQ(1, Factorial(-1));
        EXPECT_GT(Factorial(-10), 0);
        */
        // <TechnicalDetails>
        //
        // EXPECT_EQ(expected, actual) is the same as
        //
        //   EXPECT_TRUE((expected) == (actual))
        //
        // except that it will print both the expected value and the actual
        // value when the assertion fails.  This is very helpful for
        // debugging.  Therefore in this case EXPECT_EQ is preferred.
        //
        // On the other hand, EXPECT_TRUE accepts any Boolean expression,
        // and is thus more general.
        //
        // </TechnicalDetails>
    }

// Tests factorial of 0.
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
