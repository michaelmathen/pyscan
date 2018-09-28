#include "Test_Utilities.hpp"

#include <limits.h>
#include <random>
#include <iostream>

#include "HalfSpaceScan.hpp"

#include "gtest/gtest.h"
namespace {

//    TEST(HalfplaneScan, discrepancy) {
//
//        const static int n_size = 25;
//        const static int s_size = 100;
//        const static double rho = .01;
//        auto n_pts = pyscantest::randomPoints(n_size);
//        auto m_pts = pyscantest::randomPoints(s_size);
//        auto b_pts = pyscantest::randomPoints(s_size);
//        std::vector<double> m_weight(s_size, 1.0);
//        std::vector<double> b_weight(s_size, 1.0);
//
//        auto stat = [] (double m, double b) {
//          return fabs(m - b);
//        };
//        auto d1 = pyscan::MaxHalfPlane(n_pts, m_pts, m_weight, b_pts, b_weight, stat);
//        auto d2 = pyscan::MaxHalfPlane(n_pts, m_pts, m_weight, b_pts, b_weight, stat);
//
//        EXPECT_FLOAT_EQ(std::get<1>(d2),
//          pyscan::evaluate_line(std::get<0>(d2), m_pts, m_weight,
//                                b_pts, b_weight, stat));
//
//        // EXPECT_FLOAT_EQ(std::get<1>(d1), 0.0);
//        EXPECT_FLOAT_EQ(std::get<1>(d1),
//          pyscan::evaluate_line(std::get<0>(d1), m_pts, m_weight,
//                                b_pts, b_weight, stat));
//
//        EXPECT_FLOAT_EQ(std::get<1>(d1), std::get<1>(d2));
//
//
//        //std::cout << std::get<1>(d1) << std::endl;
//        //std::cout << std::get<0>(d1) << " " << std::get<0>(d2) << std::endl;
//        // auto line1 = std::get<0>(d1);
//        // auto line2 = std::get<0>(d2);
//        // EXPECT_TRUE(line1.approx_eq(line2));
//
//    }
//
//
//      TEST(HalfplaneScanLabels, discrepancy) {
//
//        const static int n_size = 25;
//        const static int s_size = 100;
//        const static double rho = .01;
//        auto n_pts = pyscantest::randomPoints(n_size);
//        auto m_pts = pyscantest::randomPoints(s_size);
//        auto b_pts = pyscantest::randomPoints(s_size);
//        auto m_labels = pyscantest::randomLabels(s_size, 30);
//        auto b_labels = pyscantest::randomLabels(s_size, 30);
//
//        std::vector<double> m_weight(s_size, 1.0);
//        std::vector<double> b_weight(s_size, 1.0);
//
//        auto stat = [] (double m, double b) {
//          return fabs(m - b);
//        };
//        pyscan::Point<> d1, d2;
//        double d1value, d2value;
//        std::tie(d1, d1value) = pyscan::max_halfplane_labeled(n_pts, m_pts, m_weight, m_labels, b_pts, b_weight, b_labels, stat);
//        std::tie(d2, d2value) = pyscan::max_halfplane_simple_labeled(n_pts, m_pts, m_weight, m_labels, b_pts, b_weight, b_labels, stat);
//
//        EXPECT_FLOAT_EQ(d2value,
//          pyscan::evaluate_line_labeled(d2, m_pts, m_weight, m_labels,
//                                b_pts, b_weight, b_labels, stat));
//
//        // EXPECT_FLOAT_EQ(std::get<1>(d1), 0.0);
//        EXPECT_FLOAT_EQ(d1value,
//          pyscan::evaluate_line_labeled(d1, m_pts, m_weight, m_labels,
//                                b_pts, b_weight, b_labels, stat));
//
//        EXPECT_FLOAT_EQ(d1value, d2value);
//
//
//        //std::cout << std::get<1>(d1) << std::endl;
//        //std::cout << std::get<0>(d1) << " " << std::get<0>(d2) << std::endl;
//        // auto line1 = std::get<0>(d1);
//        // auto line2 = std::get<0>(d2);
//        // EXPECT_TRUE(line1.approx_eq(line2));
//
//    }
}