//
// Created by mmath on 10/2/17.
//


//#include "../src/RectangleScan.hpp"
#include "DiskScan.hpp"
#include "Statistics.hpp"
#include "Test_Utilities.hpp"

#include <tuple>
#include <limits>
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


    static auto scan = [](double m, double m_total, double b, double b_total) {
        return fabs(m / m_total - b / b_total);
    };

    TEST(DiskTest, matching) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomWPoints(s_size);
        auto b_pts = pyscantest::randomWPoints(s_size);

        double d1value, d2value;
        pyscan::Disk d1;
        pyscan::Disk d2;

        std::tie(d2, d2value) = pyscan::max_disk(n_pts, m_pts, b_pts, scan);
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        std::tie(d1, d1value) = pyscan::max_disk_simple(n_pts, m_pts, b_pts, scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));

    }


    TEST(max_disk_scale, matching) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomWPoints(s_size);
        auto b_pts = pyscantest::randomWPoints(s_size);

        double d1value;
        pyscan::Disk d1;
        std::tie(d1, d1value) = max_disk_scale(n_pts, m_pts, b_pts, 1 / 32.0, 2/32.0 , scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));


    }



    TEST(DiskScan2, matching) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomWPoints(s_size);
        auto b_pts = pyscantest::randomWPoints(s_size);

        pyscan::Disk d1, d2;
        double d1value, d2value;
        std::tie(d1, d1value) = max_disk_lift(n_pts, m_pts, b_pts, scan);


        std::tie(d2, d2value) = max_disk_simple(n_pts, m_pts, b_pts, scan);

        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d1value, d2value);

    }

    void label_test(size_t n_size, size_t s_size, size_t labels) {
      auto n_pts = pyscantest::randomPoints(n_size);
      auto m_pts = pyscantest::randomLPoints(s_size, labels);
      auto b_pts = pyscantest::randomLPoints(s_size, labels);

      pyscan::Disk d1, d2;
      double d1value, d2value;
      std::tie(d1, d1value) = pyscan::max_disk_labeled(n_pts, m_pts, b_pts, scan);
      EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));

      std::tie(d2, d2value) = pyscan::max_disk_simple_labeled(n_pts, m_pts, b_pts, scan);
      EXPECT_FLOAT_EQ(d1value, d2value);

      EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
    }

    TEST(DiskTestLabel, allsamelabel) {
      label_test(25, 100, 1);
    }

    TEST(DiskTestLabel, mixedlabels) {
      label_test(25, 100, 5);
    }

//    TEST(DiskScan2, allsamelabel) {
//        label_test2(25, 100, 1);
//    }
//
//    TEST(DiskScan2, mixedlabels) {
//        label_test2(25, 100, 5);
//    }

    TEST(DiskTestLabel, alluniquelabel) {

        const static int n_size = 25;
        const static int s_size = 100;

        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_lpts = pyscantest::randomLPointsUnique(s_size);
        auto b_lpts = pyscantest::randomLPointsUnique(s_size);

        auto m_pts = pyscantest::removeLabels(m_lpts);
        auto b_pts = pyscantest::removeLabels(b_lpts);

        pyscan::Disk d1, d2, d3, d4;
        double d1value, d2value, d3value, d4value;
        std::tie(d1, d1value) = pyscan::max_disk_labeled(n_pts, m_lpts, b_lpts, scan);
        std::tie(d2, d2value) = pyscan::max_disk(n_pts, m_pts, b_pts, scan);

        std::tie(d3, d3value) = pyscan::max_disk_simple_labeled(n_pts, m_lpts, b_lpts, scan);
        std::tie(d4, d4value) = pyscan::max_disk_simple(n_pts, m_pts, b_pts, scan);


        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d1value, d2value);
        EXPECT_FLOAT_EQ(d1value, d3value);
        EXPECT_FLOAT_EQ(d2value, d4value);
        EXPECT_FLOAT_EQ(d3value, d4value);

    }


    TEST(DiskTestLabelSlow, matching) {

        const static int n_size = 25;
        const static int s_size = 100;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_lpts = pyscantest::randomLPointsUnique(s_size);
        auto b_lpts = pyscantest::randomLPointsUnique(s_size);

        auto m_pts = pyscantest::removeLabels(m_lpts);
        auto b_pts = pyscantest::removeLabels(b_lpts);

        pyscan::Disk d1, d2;
        double d1value, d2value;
        std::tie(d1, d1value) = pyscan::max_disk_simple_labeled(n_pts, m_lpts, b_lpts, scan);
        std::tie(d2, d2value) = pyscan::max_disk_simple(n_pts, m_pts, b_pts, scan);

        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d1value, d2value);

    }


}
