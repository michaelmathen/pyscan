
#include "Test_Utilities.hpp"

#include <limits.h>
#include <random>
#include <iostream>

#include "Range.hpp"
#include "HalfSpaceScan.hpp"

#include "gtest/gtest.h"
namespace {


    pyscan::discrepancy_func_t stat = [](double m, double m_total, double b, double b_total) {
        return fabs(m / m_total - b / b_total);
    };

    TEST(max_halfplane, discrepancy) {

        const static int n_size = 25;
        const static int s_size = 100;
        auto n_pts = pyscantest::randomPoints2(n_size);
        auto m_pts = pyscantest::randomWPoints2(s_size);
        auto b_pts = pyscantest::randomWPoints2(s_size);

        //std::cout << m_pts << std::endl;
        auto d1 = pyscan::max_halfplane(n_pts, m_pts, b_pts, stat);
        auto d2 = pyscan::max_halfplane_simple(n_pts, m_pts, b_pts, stat);

        EXPECT_FLOAT_EQ(std::get<1>(d2), pyscan::evaluate_range(std::get<0>(d2), m_pts, b_pts, stat));

        EXPECT_FLOAT_EQ(std::get<1>(d1),
          pyscan::evaluate_range(std::get<0>(d1), m_pts, b_pts, stat));

        EXPECT_FLOAT_EQ(std::get<1>(d1), std::get<1>(d2));
    }

    TEST(max_halfspace, discrepancy) {

        const static int n_size = 25;
        const static int s_size = 100;
        auto n_pts = pyscantest::randomPoints3(n_size);
        auto m_pts = pyscantest::randomWPoints3(s_size);
        auto b_pts = pyscantest::randomWPoints3(s_size);

        auto d1 = pyscan::max_halfspace(n_pts, m_pts, b_pts, stat);
        auto d2 = pyscan::max_halfspace_simple(n_pts, m_pts, b_pts, stat);

        EXPECT_FLOAT_EQ(std::get<1>(d2), pyscan::evaluate_range(std::get<0>(d2), m_pts, b_pts, stat));

        EXPECT_FLOAT_EQ(std::get<1>(d1), pyscan::evaluate_range(std::get<0>(d1), m_pts, b_pts, stat));

        EXPECT_FLOAT_EQ(std::get<1>(d1), std::get<1>(d2));
    }


      TEST(max_halfplane_labeled, discrepancy) {

        const static int n_size = 25;
        const static int s_size = 100;
        auto n_pts = pyscantest::randomPoints2(n_size);
        auto m_pts = pyscantest::randomLPoints2(s_size, 10);
        auto b_pts = pyscantest::randomLPoints2(s_size, 10);


        auto [d1, d1value] = pyscan::max_halfplane_labeled(n_pts, m_pts, b_pts, stat);
        auto [d2, d2value] = pyscan::max_halfplane_labeled_simple(n_pts, m_pts, b_pts, stat);

        EXPECT_FLOAT_EQ(d2value, pyscan::evaluate_range(d2, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(d1value, pyscan::evaluate_range(d1, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(d1value, d2value);

    }

    TEST(max_halfspace_labeled, discrepancy) {

        const static int n_size = 25;
        const static int s_size = 100;
        auto n_pts = pyscantest::randomPoints3(n_size);
        auto m_pts = pyscantest::randomLPoints3(s_size, 10);
        auto b_pts = pyscantest::randomLPoints3(s_size, 10);


        auto [d1, d1value] = pyscan::max_halfspace_labeled(n_pts, m_pts, b_pts, stat);
        auto [d2, d2value] = pyscan::max_halfspace_labeled_simple(n_pts, m_pts, b_pts, stat);

        EXPECT_FLOAT_EQ(d2value, pyscan::evaluate_range(d2, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(d1value, pyscan::evaluate_range(d1, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(d1value, d2value);

    }

}