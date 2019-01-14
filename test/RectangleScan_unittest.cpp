

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

    static auto scan = [](double m, double m_total, double b, double b_total) {
        return fabs(m / m_total - b / b_total);
    };

    void rect_label_test(size_t n_size, size_t s_size, size_t labels) {
        auto m_pts = pyscantest::randomLPoints2(s_size, labels);
        auto b_pts = pyscantest::randomLPoints2(s_size, labels);
        //std::cout << m_pts << std::endl;
        auto[d1, d1value] = pyscan::max_rect_labeled(n_size, std::numeric_limits<double>::infinity(), m_pts, b_pts, scan);
        std::cout << d1.upX() << " " << d1.upY() << " " << d1.lowX() << " " << d1.lowY() << std::endl;
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        //auto[d2, d2value] = pyscan::max_range3<pyscan::Rectangle, pyscan::LPoint>(n_pts, m_pts, b_pts, scan);

        //std::cout << d2value << " " << d1value << std::endl;
        //EXPECT_FLOAT_EQ(d1value, d2value);

        //EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
    }

    TEST(max_rectangle_labels, matching) {
        rect_label_test(2, 2000, 100);
    }

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
