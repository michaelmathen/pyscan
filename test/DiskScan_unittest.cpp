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



    TEST(DiskTest, matching) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomWPoints(s_size);
        auto b_pts = pyscantest::randomWPoints(s_size);
        auto scan = [](double m, double b) {
            return fabs(m - b);
        };

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

        double m_total = pyscan::computeTotal(m_pts.begin(), m_pts.end());
        double b_total = pyscan::computeTotal(b_pts.begin(), b_pts.end());
        auto scan = [](double m, double b) {
            return fabs(m - b);
        };

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
        auto scan = [](double m, double b) {
            return fabs(m - b);
        };

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
        auto scan = [](double m, double b) {
            return fabs(m - b);
        };

        double d1value;
        pyscan::Disk d1;
        std::tie(d1, d1value) = disk_scan_scale(n_pts, m_pts, b_pts, 32, scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));


    }
//
//
//    TEST(DiskScanCachedLabels, matching) {
//
//        const static int n_size = 50;
//        const static int s_size = 1000;
//        auto n_pts = pyscantest::randomPoints(n_size);
//        auto m_pts = pyscantest::randomLPoints(s_size, 50);
//        auto b_pts = pyscantest::randomLPoints(s_size, 50);
//        auto scan = [](double m, double b) {
//            return fabs(m - b);
//        };
//
//        double d1value;
//        pyscan::Disk d1;
//        //pyscan::Disk d2;
//        //std::tie(d2, d2value) = pyscan::diskScan(n_pts, m_pts, b_pts, scan);
//        //EXPECT_FLOAT_EQ(d2value, evaluate_range(m_pts, b_pts, d2, scan));
//        //std::tie(d1, d1value) = pyscan::disk_scan_simp(n_pts, m_pts, b_pts, scan);
//
//        std::tie(d1, d1value) = disk_scan_scale(n_pts, m_pts,
//                                                b_pts,
//                                                32, scan);
//        EXPECT_FLOAT_EQ(d1value, evaluate_range(m_pts, b_pts, d1, scan));
//
////        EXPECT_FLOAT_EQ(d1.getA(), d2.getA());
////        EXPECT_FLOAT_EQ(d1.getB(), d2.getB());
////        EXPECT_FLOAT_EQ(d1.getR(), d2.getR());
//
//    }
//
//    TEST(DiskScan2, matching) {
//
//        const static int n_size = 50;
//        const static int s_size = 1000;
//        auto n_pts = pyscantest::randomPoints(n_size);
//        auto m_pts = pyscantest::randomWPoints(s_size);
//        auto b_pts = pyscantest::randomWPoints(s_size);
//
//
//        auto stat = [] (double m, double b) {
//            return fabs(m - b);
//        };
//
//        pyscan::Disk d1, d2;
//        double d1value, d2value;
//        std::tie(d1, d1value) = pyscan::max_disk(n_pts, m_pts, b_pts, stat);
//
//
//        std::tie(d2, d2value) = pyscan::diskScanSlow(n_pts, m_pts, b_pts, stat);
//
//        EXPECT_FLOAT_EQ(d2value,
//                        pyscan::evaluate_range(m_pts, b_pts, d2, stat));
//
//        // EXPECT_FLOAT_EQ(std::get<1>(d1), 0.0);
//        EXPECT_FLOAT_EQ(d1value,
//                        pyscan::evaluate_range(m_pts,
//                                              b_pts, d1, stat));
//
//        EXPECT_FLOAT_EQ(d1value, d2value);
//
//    }
//
//    //updateCount tests.
//    TEST(updateCounts, adding) {
//        std::unordered_map<size_t, size_t> curr_counts;
//        curr_counts[0] = 2;
//        curr_counts[1] = 1;
//        curr_counts[3] = 1;
//
//        pyscan::crescent_t adding_set;
//        pyscan::crescent_t removing_set;
//        adding_set.emplace_back(0, .1); // shouldn't change value.
//        adding_set.emplace_back(1, .3); // shouldn't change the value.
//        adding_set.emplace_back(4, .2); // should increase the value.
//        double result_v = pyscan::updateCounts(curr_counts, adding_set, removing_set);
//
//        EXPECT_FLOAT_EQ(.2, result_v);
//        EXPECT_EQ(curr_counts[0], 3);
//        EXPECT_EQ(curr_counts[1], 2);
//        EXPECT_EQ(curr_counts[3], 1);
//        EXPECT_EQ(curr_counts[4], 1);
//    }
//
//    TEST(updateCounts, removing) {
//        std::unordered_map<size_t, size_t> curr_counts;
//        curr_counts[0] = 3;
//        curr_counts[1] = 2;
//        curr_counts[3] = 1;
//        curr_counts[4] = 1;
//
//        pyscan::crescent_t adding_set;
//        pyscan::crescent_t removing_set;
//        removing_set.emplace_back(0, .1); // shouldn't remove this element.
//        removing_set.emplace_back(1, .3); // shouldn't remove this element.
//        removing_set.emplace_back(3, .3); //should remove this element.
//        double result_v = pyscan::updateCounts(curr_counts, adding_set, removing_set);
//
//        EXPECT_FLOAT_EQ(-.3, result_v);
//        EXPECT_EQ(curr_counts[0], 2);
//        EXPECT_EQ(curr_counts[1], 1);
//        EXPECT_EQ(curr_counts.find(3), curr_counts.end());
//        EXPECT_EQ(curr_counts[4], 1);
//    }
//
//    TEST(updateCounts, addremove) {
//            std::unordered_map<size_t, size_t> curr_counts;
//            curr_counts[0] = 3;
//            curr_counts[1] = 2;
//            curr_counts[3] = 1;
//
//            pyscan::crescent_t adding_set;
//            pyscan::crescent_t removing_set;
//            adding_set.emplace_back(4, .2); // should increase the value
//            removing_set.emplace_back(4, .2);
//            double result_v = pyscan::updateCounts(curr_counts, adding_set, removing_set);
//
//            EXPECT_FLOAT_EQ(0, result_v);
//            EXPECT_EQ(curr_counts[0], 3);
//            EXPECT_EQ(curr_counts[1], 2);
//            EXPECT_EQ(curr_counts[3], 1);
//            EXPECT_EQ(curr_counts.find(4), curr_counts.end());
//    }
//
//    TEST(updateCounts, removeadd) {
//            std::unordered_map<size_t, size_t> curr_counts;
//            curr_counts[0] = 3;
//            curr_counts[1] = 2;
//            curr_counts[3] = 1;
//            curr_counts[4] = 1;
//
//            pyscan::crescent_t adding_set;
//            pyscan::crescent_t removing_set;
//            adding_set.emplace_back(4, .2); // should increase the value
//            removing_set.emplace_back(4, .2);
//            double result_v = pyscan::updateCounts(curr_counts, adding_set, removing_set);
//
//            EXPECT_FLOAT_EQ(0, result_v);
//            EXPECT_EQ(curr_counts[0], 3);
//            EXPECT_EQ(curr_counts[1], 2);
//            EXPECT_EQ(curr_counts[3], 1);
//            EXPECT_EQ(curr_counts[4], 1);
//    }


//    void label_test2(size_t n_size, size_t s_size, size_t labels) {
//        const static double rho = .01;
//        auto n_pts = pyscantest::randomPoints(n_size);
//        auto m_pts = pyscantest::randomPoints(s_size);
//        auto b_pts = pyscantest::randomPoints(s_size);
//
//        std::vector<double> m_weight(s_size, 1.0);
//        std::vector<double> b_weight(s_size, 1.0);
//        auto m_labels = pyscantest::randomLabels(s_size, 30);
//        auto b_labels = pyscantest::randomLabels(s_size, 30);
//        auto m_lpts = pyscantest::addLabels(pyscantest::addWeights(m_pts), m_labels);
//        auto b_lpts = pyscantest::addLabels(pyscantest::addWeights(b_pts), b_labels);
//        auto scan = [](double m, double b) {
//            return fabs(m - b);
//        };
//
//        pyscan::Disk d1, d2;
//        double d1value, d2value;
//        std::tie(d1, d1value) = pyscan::diskScanSlowLabels(n_pts, m_lpts, b_lpts, scan);
//
//        EXPECT_FLOAT_EQ(d1value, evaluate_range(m_lpts, b_lpts, d1, scan));
//        std::tie(d2, d2value) = pyscan::max_disk_labeled(n_pts, m_pts, m_weight, m_labels, b_pts,b_weight, b_labels, scan);
//        EXPECT_FLOAT_EQ(d1value, d2value);
//        EXPECT_FLOAT_EQ(d2value, evaluate_range(m_lpts, b_lpts, d2, scan));
//        EXPECT_FLOAT_EQ(d1.getA(), d2.getA());
//        EXPECT_FLOAT_EQ(d1.getB(), d2.getB());
//        EXPECT_FLOAT_EQ(d1.getR(), d2.getR());
//    }

//    void label_test(size_t n_size, size_t s_size, size_t labels) {
//      auto n_pts = pyscantest::randomPoints(n_size);
//      auto m_pts = pyscantest::randomLPoints(s_size, labels);
//      auto b_pts = pyscantest::randomLPoints(s_size, labels);
//      auto scan = [](double m, double b) {
//          return fabs(m - b);
//      };
//
//      pyscan::Disk d1, d2;
//      double d1value, d2value;
//      std::tie(d1, d1value) = pyscan::diskScanSlowLabels(n_pts, m_pts, b_pts, scan);
//
//      EXPECT_FLOAT_EQ(d1value, evaluate_range(m_pts, b_pts, d1, scan));
//      std::tie(d2, d2value) = pyscan::diskScanLabels(n_pts, m_pts, b_pts, scan);
//      EXPECT_FLOAT_EQ(d1value, d2value);
//      EXPECT_FLOAT_EQ(d2value, evaluate_range(m_pts, b_pts, d2, scan));
//      EXPECT_FLOAT_EQ(d1.getA(), d2.getA());
//      EXPECT_FLOAT_EQ(d1.getB(), d2.getB());
//      EXPECT_FLOAT_EQ(d1.getR(), d2.getR());
//    }
//
//    TEST(DiskTestLabel, allsamelabel) {
//      label_test(25, 100, 1);
//    }
//
//    TEST(DiskTestLabel, mixedlabels) {
//      label_test(25, 100, 5);
//    }
//
////    TEST(DiskScan2, allsamelabel) {
////        label_test2(25, 100, 1);
////    }
////
////    TEST(DiskScan2, mixedlabels) {
////        label_test2(25, 100, 5);
////    }
//
//    TEST(DiskTestLabel, alluniquelabel) {
//
//        const static int n_size = 25;
//        const static int s_size = 100;
//
//        auto f = [&](double m, double b) {
//            return fabs(m - b);
//        };
//        auto n_pts = pyscantest::randomPoints(n_size);
//        auto m_lpts = pyscantest::randomLPointsUnique(s_size);
//        auto b_lpts = pyscantest::randomLPointsUnique(s_size);
//
//        auto m_pts = pyscantest::removeLabels(m_lpts);
//        auto b_pts = pyscantest::removeLabels(b_lpts);
//
//        pyscan::Disk d1, d2, d3, d4;
//        double d1value, d2value, d3value, d4value;
//        std::tie(d1, d1value) = pyscan::diskScanLabels(n_pts, m_lpts, b_lpts, f);
//        std::tie(d2, d2value) = pyscan::diskScan(n_pts, m_pts, b_pts, f);
//
//        std::tie(d3, d3value) = pyscan::diskScanSlowLabels(n_pts, m_lpts, b_lpts, f);
//        std::tie(d4, d4value) = pyscan::diskScanSlow(n_pts, m_pts, b_pts, f);
//
//
//        EXPECT_FLOAT_EQ(d1value, evaluate_range(m_pts, b_pts, d1, f));
//        EXPECT_FLOAT_EQ(d2value, evaluate_range(m_pts, b_pts, d2, f));
//        EXPECT_FLOAT_EQ(d1value, d2value);
//        EXPECT_FLOAT_EQ(d1value, d3value);
//        EXPECT_FLOAT_EQ(d2value, d4value);
//        EXPECT_FLOAT_EQ(d3value, d4value);
//
//        EXPECT_FLOAT_EQ(d1.getA(), d2.getA());
//        EXPECT_FLOAT_EQ(d1.getB(), d2.getB());
//        EXPECT_FLOAT_EQ(d1.getR(), d2.getR());
//    }
//
//
//    TEST(DiskTestLabelSlow, matching) {
//
//        const static int n_size = 25;
//        const static int s_size = 100;
//        auto n_pts = pyscantest::randomPoints(n_size);
//        auto m_lpts = pyscantest::randomLPointsUnique(s_size);
//        auto b_lpts = pyscantest::randomLPointsUnique(s_size);
//
//        auto m_pts = pyscantest::removeLabels(m_lpts);
//        auto b_pts = pyscantest::removeLabels(b_lpts);
//
//        auto f = [&](double m, double b) {
//            return fabs(m - b);
//        };
//
//        pyscan::Disk d1, d2;
//        double d1value, d2value;
//        std::tie(d1, d1value) = pyscan::diskScanSlowLabels(n_pts, m_lpts, b_lpts, f);
//        std::tie(d2, d2value) = pyscan::diskScanSlow(n_pts, m_pts, b_pts, f);
//
//        EXPECT_FLOAT_EQ(d1value, evaluate_range(m_pts, b_pts, d1, f));
//        EXPECT_FLOAT_EQ(d2value, evaluate_range(m_pts, b_pts, d2, f));
//        EXPECT_FLOAT_EQ(d1value, d2value);
//        EXPECT_FLOAT_EQ(d1.getA(), d2.getA());
//        EXPECT_FLOAT_EQ(d1.getB(), d2.getB());
//        EXPECT_FLOAT_EQ(d1.getR(), d2.getR());
//    }
//
//    TEST(computeLabelTotalF, inserting) {
//      std::unordered_map<size_t, size_t> curr_counts;
//      curr_counts[2] = 1;
//
//      std::vector<pyscan::LPoint<>> bunch_of_points = {
//        pyscan::LPoint<>(0, .1, 0.0, .1, 1.0),
//        pyscan::LPoint<>(0, .1, .1, .2, 1.0),
//        pyscan::LPoint<>(0, .1, .1, .2, 1.0),
//        pyscan::LPoint<>(1, .1, .1, .2, 1.0),
//        pyscan::LPoint<>(2, .1, .1, .2, 1.0),
//      };
//
//      double fv = pyscan::computeLabelTotal(bunch_of_points.begin(), bunch_of_points.end(),
//                                  curr_counts);
//
//      EXPECT_FLOAT_EQ(fv, .2);
//      EXPECT_EQ(curr_counts[0], 3);
//      EXPECT_EQ(curr_counts[1], 1);
//      EXPECT_EQ(curr_counts[2], 2);
//
//      curr_counts = std::unordered_map<size_t, size_t>();
//      fv = pyscan::computeLabelTotalF(bunch_of_points.begin(), bunch_of_points.end(),
//                                      curr_counts,
//                                      [&] (pyscan::Point<> const& pt) {
//                                        return pt(0) >.05;
//                                      });
//      EXPECT_FLOAT_EQ(fv, .3);
//      EXPECT_EQ(curr_counts[0], 2);
//      EXPECT_EQ(curr_counts[1], 1);
//      EXPECT_EQ(curr_counts[2], 1);
//    }

}
