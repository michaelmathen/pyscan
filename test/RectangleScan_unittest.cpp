


#include "RectangleScan.hpp"
#include "Utilities.hpp"

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


    std::tuple<std::vector<size_t>, pyscan::epoint_list_t, pyscan::epoint_list_t> initialize_pts() {
        const static int n_size = 50;
        const static int s_size = 1000;
        auto m_wpts = pyscantest::randomWPoints2(s_size);
        auto b_wpts = pyscantest::randomWPoints2(s_size);

        auto [m_pts, b_pts] = pyscan::to_epoints(m_wpts, b_wpts);

        auto n_pts = pyscantest::samplewr(b_pts, n_size);

        std::vector<size_t> divisions;
        for (auto& p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());
        return std::make_tuple(divisions, m_pts, b_pts);
    }

    TEST(SlabTree, init) {

        auto [divisions, m_pts, b_pts] = initialize_pts();

        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);
        auto root = tree.get_root();

        std::vector<decltype(root)> curr_stack;
        curr_stack.push_back(root);
        std::vector<decltype(root)> leaves;
        while (!curr_stack.empty()) {
            auto el = curr_stack.back();

            ASSERT_LE(el->bottom_y, el->top_y);

            curr_stack.pop_back();
            if (el->up != nullptr) {
                ASSERT_EQ(el->top_y, el->up->top_y);
                ASSERT_LE(el->up->bottom_y, el->top_y);
                ASSERT_LE(el->bottom_y, el->up->bottom_y);
                curr_stack.push_back(el->up);
            }
            if (el->down != nullptr) {
                ASSERT_EQ(el->bottom_y, el->down->bottom_y);
                ASSERT_LE(el->down->top_y, el->top_y);
                ASSERT_LE(el->bottom_y, el->down->top_y);
                curr_stack.push_back(el->down);
            }
            if (el->has_mid()) {
                ASSERT_EQ(el->down->top_y, el->up->bottom_y);
            }
            if (!el->up  && !el->down){
                leaves.push_back(el);
            }
        }

        std::vector<double> divisions2;
        for (auto& leaf : leaves) {
            divisions2.push_back(leaf->bottom_y);
        }

        using pyscan::operator<<;
        divisions2.push_back(leaves.back()->top_y);
        ASSERT_EQ(divisions.size(), divisions2.size());
        for (size_t i = 0; i < divisions.size(); i++) {
            ASSERT_EQ(divisions2[i], divisions[i]);
        }

        size_t curr_i = 0;
        for (auto& leaf : leaves) {
            ASSERT_EQ(divisions[curr_i], leaf->bottom_y);
            ASSERT_EQ(divisions[curr_i + 1], leaf->top_y);
            curr_i++;
        }

    }


    TEST(Slab, get_contains) {
        auto [divisions, m_pts, b_pts] = initialize_pts();

        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);
        size_t x_min = std::min(m_pts[0](0), m_pts[1](0));
        size_t x_max = std::max(m_pts[0](0), m_pts[1](0));

        pyscan::ERectangle rect(x_max, divisions[5], x_min, divisions[2]);
        auto slab_containing = tree.get_containing(rect);
        ASSERT_GT(slab_containing->top_y, rect.upY());
        ASSERT_LE(slab_containing->bottom_y, rect.lowY());

        ASSERT_GE(slab_containing->get_mid(), rect.lowY());
        ASSERT_LT(slab_containing->get_mid(), rect.upY());

    }

    TEST(SlabTree, measure_interval) {
        auto [divisions, m_pts, b_pts] = initialize_pts();

        pyscan::SlabTree tree(divisions, m_pts, b_pts, true, 1.0);
        auto root = tree.get_root();
        if (!root || !root->down) {
            return;
        }
        double val = root->measure_interval(20, 4, 1.0, 1.0);
        pyscan::ERectangle rect1(.6, root->top_y, .2, root->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect1, m_pts) + pyscan::range_weight(rect1, b_pts));

        val = root->down->measure_interval(20, 4, 1.0, 1.0);
        ASSERT_GE(root->down->top_y, root->down->bottom_y);
        pyscan::ERectangle rect2(.6, root->down->top_y, .2, root->down->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect2, m_pts) + pyscan::range_weight(rect2, b_pts));

        val = root->up->measure_interval(20, 4, 1.0, 1.0);
        ASSERT_GE(root->up->top_y, root->up->bottom_y);
        pyscan::ERectangle rect3(.6, root->up->top_y, .2, root->up->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect3, m_pts) + pyscan::range_weight(rect3, b_pts));


        while (root->up != nullptr) {
            root = root->up;
        }
        val = root->measure_interval(20, 4, 1.0, 1.0);
        ASSERT_GE(root->top_y, root->bottom_y);
        pyscan::ERectangle rect4(.6, root->top_y, .2, root->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect4, m_pts) + pyscan::range_weight(rect4, b_pts));

    }

    TEST(SlabTree, measure_rect) {

        const static int n_size = 20;
        const static int s_size = 1000;
        auto m_wpts = pyscantest::randomWPoints2(s_size);
        auto b_wpts = pyscantest::randomWPoints2(s_size);

        auto [m_pts, b_pts] = pyscan::to_epoints(m_wpts, b_wpts);

        auto n_pts = pyscantest::samplewr(b_pts, n_size);

        std::vector<size_t> divisions;
        for (auto& p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());

        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);
        for (auto it1 = n_pts.begin(); it1 != n_pts.end() - 3; ++it1) {
            for (auto it2 = it1 + 1; it2 != n_pts.end() - 3; ++it2) {
                for (auto it3 = it2 + 1; it3 != n_pts.end() - 2; ++it3) {
                    for (auto it4 = it3 + 1; it4 != n_pts.end() - 1; ++it4) {
                        pyscan::ERectangle rect(*it1, *it2, *it3, *it4);
                        double val = tree.measure_rect(rect, 1.0, 1.0);
                        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect, m_pts) + pyscan::range_weight(rect, b_pts));
                    }
                }
            }
        }
    }

    TEST(SlabTree, split_offsets) {
        /*
         * Test that the split offsets were created correctly in the SlabTree.
         */
        const static int n_size = 20;
        const static int s_size = 1000;
        auto m_wpts = pyscantest::randomWPoints2(s_size);
        auto b_wpts = pyscantest::randomWPoints2(s_size);

        auto [m_pts, b_pts] = pyscan::to_epoints(m_wpts, b_wpts);

        auto n_pts = pyscantest::samplewr(b_pts, n_size);

        std::vector<size_t> divisions;
        for (auto& p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());

        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);
        for (auto it1 = n_pts.begin(); it1 != n_pts.end() - 3; ++it1) {
            for (auto it2 = it1 + 1; it2 != n_pts.end() - 3; ++it2) {
                for (auto it3 = it2 + 1; it3 != n_pts.end() - 2; ++it3) {
                    for (auto it4 = it3 + 1; it4 != n_pts.end() - 1; ++it4) {
                        pyscan::ERectangle rect(*it1, *it2, *it3, *it4);
                        double val = tree.measure_rect(rect, 1.0, 1.0);
                        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect, m_pts) + pyscan::range_weight(rect, b_pts));
                    }
                }
            }
        }
    }

    TEST(SlabTree, max_rectangle) {

        const static int n_size = 20;
        const static int s_size = 1000;
        auto m_wpts = pyscantest::randomWPoints2(s_size);
        auto b_wpts = pyscantest::randomWPoints2(s_size);

        auto [m_pts, b_pts] = pyscan::to_epoints(m_wpts, b_wpts);

        auto n_pts = pyscantest::samplewr(b_pts, n_size);

        std::vector<size_t> divisions;
        for (auto& p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());

        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);

        auto [rect, val] = tree.max_rectangle(1.0, -1.0);
        std::cout << val << std::endl;
        std::cout << rect.toString() << std::endl;
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect, m_pts) - pyscan::range_weight(rect, b_pts));


    }
}
