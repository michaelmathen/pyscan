


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

    TEST(Slab, measure_interval) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints2(n_size);
        auto m_pts = pyscantest::randomWPoints2(s_size);
        auto b_pts = pyscantest::randomWPoints2(s_size);
        std::vector<double> divisions;
        for (auto& p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());

        pyscan::SlabTree tree(divisions, m_pts, b_pts, 1.0);
        auto root = tree.get_root();
        if (!root || !root->down) {
            return;
        }
        double val = root->measure_interval(.6, .2, 1.0, 1.0);
        pyscan::Rectangle rect1(.6, root->top_y, .2, root->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect1, m_pts) + pyscan::range_weight(rect1, b_pts));

        val = root->down->measure_interval(.6, .2, 1.0, 1.0);
        ASSERT_GE(root->down->top_y, root->down->bottom_y);
        pyscan::Rectangle rect2(.6, root->down->top_y, .2, root->down->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect2, m_pts) + pyscan::range_weight(rect2, b_pts));

        val = root->up->measure_interval(.6, .2, 1.0, 1.0);
        ASSERT_GE(root->up->top_y, root->up->bottom_y);
        pyscan::Rectangle rect3(.6, root->up->top_y, .2, root->up->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect3, m_pts) + pyscan::range_weight(rect3, b_pts));
    }

    TEST(SlabTree, approx) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints2(n_size);
        auto m_pts = pyscantest::randomWPoints2(s_size);
        auto b_pts = pyscantest::randomWPoints2(s_size);
        std::vector<double> divisions;
        for (auto& p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());

        pyscan::SlabTree tree(divisions, m_pts, b_pts, 1.0);
        for (auto it1 = n_pts.begin(); it1 != n_pts.end() - 3; ++it1) {
            for (auto it2 = it1 + 1; it2 != n_pts.end() - 3; ++it2) {
                for (auto it3 = it2 + 1; it3 != n_pts.end() - 2; ++it3) {
                    for (auto it4 = it3 + 1; it4 != n_pts.end() - 1; ++it4) {
                        pyscan::Rectangle rect(*it1, *it2, *it3, *it4);
                        double val = tree.measure_rect(rect, 1.0, 1.0);
                        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect, m_pts) + pyscan::range_weight(rect, b_pts));
                    }
                }
            }
        }



    }


}
