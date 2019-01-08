

#include "RectangleScan.hpp"
#include "Test_Utilities.hpp"


#include "gtest/gtest.h"

namespace {

//    TEST(grid, rectangle) {
//
//        const static int s_size = 1000;
//        auto m_pts = pyscantest::randomWPoints2(s_size);
//        auto b_pts = pyscantest::randomWPoints2(s_size);
//        pyscan::Grid grid(50, m_pts, b_pts);
//
//        for (size_t i = 0; i < grid.size() - 1; i++ {
//            for (size_t j = i + 1; i < grid.size(); i++ {
//                for (size_t k = 0; k < grid.size() - 1; k++ {
//                    for (size_t l = k + 1; k < grid.size(); k++ {
//                        grid.xCoord(i)
//                    }
//                }
//            }
//        }
//    };

    TEST(max_rectangle, matching) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints2(n_size);
        auto m_pts = pyscantest::randomWPoints2(s_size);
        auto b_pts = pyscantest::randomWPoints2(s_size);

        auto scan = [](double m, double m_total, double b, double b_total) {
            return m / m_total - b / b_total;
        };
        pyscan::Grid grid(50, m_pts, b_pts);

        auto subgrid_lin = pyscan::max_subgrid_linear(grid, 1.0, -1.0);
        auto subgrid = pyscan::max_subgrid(grid, scan);

        EXPECT_FLOAT_EQ(subgrid.fValue(), subgrid_lin.fValue());
        std::cout << subgrid.fValue() << std::endl;

    }

}
