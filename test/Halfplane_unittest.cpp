#include "Test_Utilities.hpp"

#include <limits.h>
#include <random>
#include <iostream>

#include "Range.hpp"
#include "HalfSpaceScan.hpp"

#include "gtest/gtest.h"
namespace {

    TEST(HalfplaneScan, discrepancy) {

        const static int n_size = 25;
        const static int s_size = 100;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomWPoints(s_size);
        auto b_pts = pyscantest::randomWPoints(s_size);

        auto stat = [] (double m, double b) {
          return fabs(m - b);
        };
        auto d1 = pyscan::max_halfplane(n_pts, m_pts, b_pts, stat);
        auto d2 = pyscan::max_halfplane_simple(n_pts, m_pts, b_pts, stat);

        EXPECT_FLOAT_EQ(std::get<1>(d2), pyscan::evaluate_range(std::get<0>(d2), m_pts, b_pts, stat));

        EXPECT_FLOAT_EQ(std::get<1>(d1),
          pyscan::evaluate_range(std::get<0>(d1), m_pts, b_pts, stat));

        EXPECT_FLOAT_EQ(std::get<1>(d1), std::get<1>(d2));
    }


      TEST(HalfplaneScanLabels, discrepancy) {

        const static int n_size = 25;
        const static int s_size = 100;
        auto n_pts = pyscantest::randomPoints(n_size);
        auto m_pts = pyscantest::randomLPoints(s_size, 10);
        auto b_pts = pyscantest::randomLPoints(s_size, 10);


        auto stat = [] (double m, double b) {
          return fabs(m - b);
        };
        pyscan::halfspace2_t d1, d2;
        double d1value, d2value;
        std::tie(d1, d1value) = pyscan::max_halfplane_labeled(n_pts, m_pts, b_pts, stat);
        std::tie(d2, d2value) = pyscan::max_halfplane_labeled_simple(n_pts, m_pts, b_pts, stat);

        EXPECT_FLOAT_EQ(d2value, pyscan::evaluate_range(d2, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(d1value, pyscan::evaluate_range(d1, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(d1value, d2value);

    }
}