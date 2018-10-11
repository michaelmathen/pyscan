//
// Created by mmath on 10/2/17.
//


//#include "../src/RectangleScan.hpp"
#include "DiskScan.hpp"
#include "DiskScan2.hpp"
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


    static auto scan = [](double m, double b) {
        return fabs(m - b);
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

        std::tie(d2, d2value) = pyscan::disk_scan(n_pts, m_pts, b_pts, scan);
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        std::tie(d1, d1value) = pyscan::disk_scan_simple(n_pts, m_pts, b_pts, scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));

        EXPECT_FLOAT_EQ(d1.getA(), d2.getA());
        EXPECT_FLOAT_EQ(d1.getB(), d2.getB());
        EXPECT_FLOAT_EQ(d1.getR(), d2.getR());

    }

    TEST(DiskScanRestricted, matching) {

        const static int n_size = 50;
        const static int s_size = 10000;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomWPoints(s_size);
        auto b_pts = pyscantest::randomWPoints(s_size);

        double m_total = pyscan::computeTotal(m_pts);
        double b_total = pyscan::computeTotal(b_pts);

        double d2value;
        pyscan::Disk d2;
        std::tie(d2, d2value) = pyscan::disk_scan_restricted(n_pts[0], n_pts[1], n_pts, m_pts, b_pts,
                0, std::numeric_limits<double>::infinity(), m_total, b_total, scan);
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
    }

    TEST(DiskScanSimp, matching) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomWPoints(s_size);
        auto b_pts = pyscantest::randomWPoints(s_size);
        double d1value, d2value;
        pyscan::Disk d1;
        pyscan::Disk d2;
        std::tie(d2, d2value) = pyscan::disk_scan(n_pts, m_pts, b_pts, scan);
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        std::tie(d1, d1value) = pyscan::disk_scan_simple(n_pts, m_pts, b_pts, scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d1value, d2value);

    }


    TEST(DiskScanCached, matching) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomWPoints(s_size);
        auto b_pts = pyscantest::randomWPoints(s_size);

        double d1value;
        pyscan::Disk d1;
        std::tie(d1, d1value) = disk_scan_scale(n_pts, m_pts, b_pts, 32, scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));


    }


    TEST(DiskScanCachedLabels, matching) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomLPoints(s_size, 50);
        auto b_pts = pyscantest::randomLPoints(s_size, 50);

        double d1value;
        pyscan::Disk d1;
        std::tie(d1, d1value) = disk_scan_scale(n_pts, m_pts,
                                                b_pts,
                                                32, scan);
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
        std::tie(d1, d1value) = disk_scan2(n_pts, m_pts, b_pts, scan);


        std::tie(d2, d2value) = disk_scan_simple(n_pts, m_pts, b_pts, scan);

        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d1value, d2value);

    }

    //updateCount tests.
    TEST(updateCounts, adding) {
        std::unordered_map<size_t, size_t> curr_counts;
        curr_counts[0] = 2;
        curr_counts[1] = 1;
        curr_counts[3] = 1;

        pyscan::crescent_t adding_set;
        pyscan::crescent_t removing_set;
        adding_set.emplace_back(0, .1); // shouldn't change value.
        adding_set.emplace_back(1, .3); // shouldn't change the value.
        adding_set.emplace_back(4, .2); // should increase the value.
        double result_v = pyscan::updateCounts(curr_counts, adding_set, removing_set);

        EXPECT_FLOAT_EQ(.2, result_v);
        EXPECT_EQ(curr_counts[0], 3);
        EXPECT_EQ(curr_counts[1], 2);
        EXPECT_EQ(curr_counts[3], 1);
        EXPECT_EQ(curr_counts[4], 1);
    }

    TEST(updateCounts, removing) {
        std::unordered_map<size_t, size_t> curr_counts;
        curr_counts[0] = 3;
        curr_counts[1] = 2;
        curr_counts[3] = 1;
        curr_counts[4] = 1;

        pyscan::crescent_t adding_set;
        pyscan::crescent_t removing_set;
        removing_set.emplace_back(0, .1); // shouldn't remove this element.
        removing_set.emplace_back(1, .3); // shouldn't remove this element.
        removing_set.emplace_back(3, .3); //should remove this element.
        double result_v = pyscan::updateCounts(curr_counts, adding_set, removing_set);

        EXPECT_FLOAT_EQ(-.3, result_v);
        EXPECT_EQ(curr_counts[0], 2);
        EXPECT_EQ(curr_counts[1], 1);
        EXPECT_EQ(curr_counts.find(3), curr_counts.end());
        EXPECT_EQ(curr_counts[4], 1);
    }

    TEST(updateCounts, addremove) {
        std::unordered_map<size_t, size_t> curr_counts;
        curr_counts[0] = 3;
        curr_counts[1] = 2;
        curr_counts[3] = 1;

        pyscan::crescent_t adding_set;
        pyscan::crescent_t removing_set;
        adding_set.emplace_back(4, .2); // should increase the value
        removing_set.emplace_back(4, .2);
        double result_v = pyscan::updateCounts(curr_counts, adding_set, removing_set);

        EXPECT_FLOAT_EQ(0, result_v);
        EXPECT_EQ(curr_counts[0], 3);
        EXPECT_EQ(curr_counts[1], 2);
        EXPECT_EQ(curr_counts[3], 1);
        EXPECT_EQ(curr_counts.find(4), curr_counts.end());
    }

    TEST(updateCounts, removeadd) {
        std::unordered_map<size_t, size_t> curr_counts;
        curr_counts[0] = 3;
        curr_counts[1] = 2;
        curr_counts[3] = 1;
        curr_counts[4] = 1;

        pyscan::crescent_t adding_set;
        pyscan::crescent_t removing_set;
        adding_set.emplace_back(4, .2); // should increase the value
        removing_set.emplace_back(4, .2);
        double result_v = pyscan::updateCounts(curr_counts, adding_set, removing_set);

        EXPECT_FLOAT_EQ(0, result_v);
        EXPECT_EQ(curr_counts[0], 3);
        EXPECT_EQ(curr_counts[1], 2);
        EXPECT_EQ(curr_counts[3], 1);
        EXPECT_EQ(curr_counts[4], 1);
    }


    void label_test2(size_t n_size, size_t s_size, size_t labels) {
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_lpts = pyscantest::randomLPoints(s_size, labels);
        auto b_lpts = pyscantest::randomLPoints(s_size, labels);

        pyscan::Disk d1, d2;
        double d1value, d2value;
        std::tie(d1, d1value) = pyscan::disk_scan_simple_labels(n_pts, m_lpts, b_lpts, scan);

        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_lpts, b_lpts, scan));
        std::tie(d2, d2value) = pyscan::disk_scan2_labels(n_pts, m_lpts, b_lpts, scan);
        EXPECT_FLOAT_EQ(d1value, d2value);
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_lpts, b_lpts, scan));
    }

    void label_test(size_t n_size, size_t s_size, size_t labels) {
      auto n_pts = pyscantest::randomPoints(n_size);
      auto m_pts = pyscantest::randomLPoints(s_size, labels);
      auto b_pts = pyscantest::randomLPoints(s_size, labels);

      pyscan::Disk d1, d2;
      double d1value, d2value;
      std::tie(d1, d1value) = pyscan::disk_scan_simple_labels(n_pts, m_pts, b_pts, scan);
      EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));

      std::tie(d2, d2value) = pyscan::disk_scan_labels(n_pts, m_pts, b_pts, scan);
      EXPECT_FLOAT_EQ(d1value, d2value);
      std::cout << d2.getB() << " " << d2.getA() << " " << d2.getR() << std::endl;
      std::cout << d1.getB() << " " << d1.getA() << " " << d1.getR() << std::endl;

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

        auto f = [&](double m, double b) {
            return fabs(m - b);
        };
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_lpts = pyscantest::randomLPointsUnique(s_size);
        auto b_lpts = pyscantest::randomLPointsUnique(s_size);

        auto m_pts = pyscantest::removeLabels(m_lpts);
        auto b_pts = pyscantest::removeLabels(b_lpts);

        pyscan::Disk d1, d2, d3, d4;
        double d1value, d2value, d3value, d4value;
        std::tie(d1, d1value) = pyscan::disk_scan_labels(n_pts, m_lpts, b_lpts, f);
        std::tie(d2, d2value) = pyscan::disk_scan(n_pts, m_pts, b_pts, f);

        std::tie(d3, d3value) = pyscan::disk_scan_simple_labels(n_pts, m_lpts, b_lpts, f);
        std::tie(d4, d4value) = pyscan::disk_scan_simple(n_pts, m_pts, b_pts, f);


        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, f));
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, f));
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

        auto f = [&](double m, double b) {
            return fabs(m - b);
        };

        pyscan::Disk d1, d2;
        double d1value, d2value;
        std::tie(d1, d1value) = pyscan::disk_scan_simple_labels(n_pts, m_lpts, b_lpts, f);
        std::tie(d2, d2value) = pyscan::disk_scan_simple(n_pts, m_pts, b_pts, f);

        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, f));
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, f));
        EXPECT_FLOAT_EQ(d1value, d2value);

    }


}
