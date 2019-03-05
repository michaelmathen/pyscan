


#include "RectangleScan.hpp"
#include "IntervalScan.hpp"
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
        auto[d1, d1value] = pyscan::max_rect_labeled(n_size, std::numeric_limits<double>::infinity(), m_pts, b_pts,
                                                     scan);
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

        auto[m_pts, b_pts, tmp1, tmp2] = pyscan::to_epoints(m_wpts, b_wpts);
        (void) tmp1;
        (void) tmp2;

        auto n_pts = pyscantest::samplewr(b_pts, n_size);

        std::vector<size_t> divisions;
        for (auto &p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());
        return std::make_tuple(divisions, m_pts, b_pts);
    }

    TEST(SlabTree, init) {

        auto[divisions, m_pts, b_pts] = initialize_pts();

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
            if (!el->up && !el->down) {
                leaves.push_back(el);
            }
        }

        std::vector<size_t> divisions2;
        for (auto &leaf : leaves) {
            divisions2.push_back(leaf->bottom_y);
        }

        using pyscan::operator<<;
        divisions2.push_back(leaves.back()->top_y);
        ASSERT_EQ(divisions.size(), divisions2.size());
        for (size_t i = 0; i < divisions.size(); i++) {
            ASSERT_EQ(divisions2[i], divisions[i]);
        }

        size_t curr_i = 0;
        for (auto &leaf : leaves) {
            ASSERT_EQ(divisions[curr_i], leaf->bottom_y);
            ASSERT_EQ(divisions[curr_i + 1], leaf->top_y);
            curr_i++;
        }

    }


    TEST(Slab, get_contains) {
        auto[divisions, m_pts, b_pts] = initialize_pts();

        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);
        size_t x_min = std::min(m_pts[0](0), m_pts[1](0));
        size_t x_max = std::max(m_pts[0](0), m_pts[1](0));

        pyscan::ERectangle rect(x_max, divisions[5], x_min, divisions[2]);
        auto slab_containing = tree.get_containing(rect.upY(), rect.lowY());
        ASSERT_GT(slab_containing->top_y, rect.upY());
        ASSERT_LE(slab_containing->bottom_y, rect.lowY());

        ASSERT_GE(slab_containing->get_mid(), rect.lowY());
        ASSERT_LT(slab_containing->get_mid(), rect.upY());

    }

    TEST(SlabTree, measure_interval) {
        auto[divisions, m_pts, b_pts] = initialize_pts();

        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);
        auto root = tree.get_root();
        if (!root || !root->down) {
            return;
        }
        double val = root->measure_interval(20, 4, 1.0, 1.0);
        pyscan::ERectangle rect1(20, root->top_y, 4, root->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect1, m_pts) + pyscan::range_weight(rect1, b_pts));

        val = root->down->measure_interval(20, 4, 1.0, 1.0);
        ASSERT_GE(root->down->top_y, root->down->bottom_y);
        pyscan::ERectangle rect2(20, root->down->top_y, 4, root->down->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect2, m_pts) + pyscan::range_weight(rect2, b_pts));

        val = root->up->measure_interval(20, 4, 1.0, 1.0);
        ASSERT_GE(root->up->top_y, root->up->bottom_y);
        pyscan::ERectangle rect3(20, root->up->top_y, 4, root->up->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect3, m_pts) + pyscan::range_weight(rect3, b_pts));


        while (root->up != nullptr) {
            root = root->up;
        }
        val = root->measure_interval(20, 4, 1.0, 1.0);
        ASSERT_GE(root->top_y, root->bottom_y);
        pyscan::ERectangle rect4(20, root->top_y, 4, root->bottom_y);
        ASSERT_FLOAT_EQ(val, pyscan::range_weight(rect4, m_pts) + pyscan::range_weight(rect4, b_pts));

    }

    TEST(SlabTree, measure_rect) {

        const static int n_size = 20;
        const static int s_size = 1000;
        auto m_wpts = pyscantest::randomWPoints2(s_size);
        auto b_wpts = pyscantest::randomWPoints2(s_size);

        auto[m_pts, b_pts, tmp1, tmp2] = pyscan::to_epoints(m_wpts, b_wpts);
        (void) tmp1;
        (void) tmp2;
        auto n_pts = pyscantest::samplewr(b_pts, n_size);

        std::vector<size_t> divisions;
        for (auto &p : n_pts) {
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
         * TODO
         */
        const static int n_size = 20;
        const static int s_size = 10;
        auto m_wpts = pyscantest::randomWPoints2(s_size);
        auto b_wpts = pyscantest::randomWPoints2(s_size);
        using pyscan::operator<<;

        auto[m_pts, b_pts, tmp1, tmp2] = pyscan::to_epoints(m_wpts, b_wpts);
        (void) tmp1;
        (void) tmp2;

        auto n_pts = pyscantest::samplewr(b_pts, n_size);

        std::vector<size_t> divisions;
        for (auto &p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());

        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);

        auto root = tree.get_root();
        std::vector<decltype(root)> curr_stack;
        curr_stack.push_back(root);
        std::vector<decltype(root)> leaves;
        while (!curr_stack.empty()) {
            auto el = curr_stack.back();
            ASSERT_TRUE(std::is_sorted(el->global_split_offset.begin(), el->global_split_offset.end()));
            curr_stack.pop_back();
            if (el->up != nullptr) {
                curr_stack.push_back(el->up);
                ASSERT_TRUE(includes(el->global_split_offset.begin(), el->global_split_offset.end(),
                                     el->up->global_split_offset.begin(), el->up->global_split_offset.end()));
            }
            if (el->down != nullptr) {
                curr_stack.push_back(el->down);
                ASSERT_TRUE(includes(el->global_split_offset.begin(), el->global_split_offset.end(),
                                     el->down->global_split_offset.begin(), el->down->global_split_offset.end()));
            }
            if (el->up != nullptr && el->down != nullptr) {
                std::vector<size_t> intersection;
                std::set_intersection(el->up->global_split_offset.begin(), el->up->global_split_offset.end(),
                                      el->down->global_split_offset.begin(), el->down->global_split_offset.end(),
                                      std::back_inserter(intersection));
                if (intersection.size() > 0) {
                    for (size_t i = 0; i < intersection.size(); i++) {
                        std::cout << intersection[i] << std::endl;
                    }
                    EXPECT_TRUE(false);
                }
            }
            if (!el->up && !el->down) {
                leaves.push_back(el);
            }
        }
    }

    TEST(SlabTree, dyadic_approx) {
        const static int n_size = 20;
        const static int s_size = 100;
        auto m_wpts = pyscantest::randomWPoints2(s_size);
        auto b_wpts = pyscantest::randomWPoints2(s_size);

        auto[m_pts, b_pts, tmp1, tmp2] = pyscan::to_epoints(m_wpts, b_wpts);
        (void) tmp1;
        (void) tmp2;
        auto n_pts = pyscantest::samplewr(b_pts, n_size);

        std::vector<size_t> divisions;
        for (auto &p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());
        auto new_end = std::unique(divisions.begin(), divisions.end());
        divisions.erase(new_end, divisions.end());
        using pyscan::operator<<;
        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);

        for (size_t i = 0; i < divisions.size(); i++){
            for (size_t j = i; j < divisions.size(); j++) {
                auto [m1, b1] = tree.dyadic_approx(divisions[j], divisions[i]);
                pyscan::epoint_list_t region_slab;
                for (auto& p : m_pts) {
                    if (divisions[i] <= p.get_y() && p.get_y() < divisions[j]) {
                        region_slab.emplace_back(p);
                    }
                }
                for (auto& p : b_pts) {
                    if (divisions[i] <= p.get_y() && p.get_y() < divisions[j]) {
                        region_slab.emplace_back(p);
                    }
                }
                pyscan::epoint_list_t region_slab2;
                for (auto& p : m1) {
                    ASSERT_LE(divisions[i], p.get_y());
                    ASSERT_LT(p.get_y(), divisions[j]);
                    region_slab2.emplace_back(p);
                }
                for (auto& p : b1) {
                    ASSERT_LE(divisions[i], p.get_y());
                    ASSERT_LT(p.get_y(), divisions[j]);
                    region_slab2.emplace_back(p);
                }

                std::sort(region_slab.begin(), region_slab.end());
                std::sort(region_slab2.begin(), region_slab2.end());

                auto it = std::unique(region_slab.begin(), region_slab.end());
                ASSERT_EQ(it - region_slab.begin(), region_slab.size());
                if (!std::includes(region_slab.begin(), region_slab.end(), region_slab2.begin(), region_slab2.end()) ||
                    !std::includes(region_slab2.begin(), region_slab2.end(), region_slab.begin(), region_slab.end())) {
                    std::cout << region_slab << std::endl;
                    std::cout << region_slab2 << std::endl;
                }
                ASSERT_TRUE(std::includes(region_slab.begin(), region_slab.end(), region_slab2.begin(), region_slab2.end()));
                ASSERT_TRUE(std::includes(region_slab2.begin(), region_slab2.end(), region_slab.begin(), region_slab.end()));

            }
        }
    }

    TEST(SlabTree, max_rectangle) {

        const static int n_size = 20;
        const static int s_size = 100;
        auto m_wpts = pyscantest::randomWPoints2(s_size);
        auto b_wpts = pyscantest::randomWPoints2(s_size);

        auto[m_pts, b_pts, tmp1, tmp2] = pyscan::to_epoints(m_wpts, b_wpts);
        (void) tmp1;
        (void) tmp2;

        auto n_pts = pyscantest::samplewr(b_pts, n_size);

        std::vector<size_t> divisions;
        for (auto &p : n_pts) {
            divisions.emplace_back(p(1));
        }
        std::sort(divisions.begin(), divisions.end());
        auto new_end = std::unique(divisions.begin(), divisions.end());
        divisions.erase(new_end, divisions.end());
        using pyscan::operator<<;
        std::cout << divisions << std::endl;
        pyscan::SlabTree tree(divisions, m_pts, b_pts, false, 1.0);

        auto[rect, val] = tree.max_rectangle(1.0, -1.0);
        auto [rect2, rect_val] = tree.max_rectangle_slow(1.0, -1.0);
        std::cout << rect.toString() << " " <<  val <<  " "
            << pyscan::range_weight(rect, m_pts) / pyscan::computeTotal(m_pts) - pyscan::range_weight(rect, b_pts) / pyscan::computeTotal(b_pts)
            << std::endl;

        std::cout << rect2.toString() << " " << rect_val << " "
            <<  pyscan::range_weight(rect2, m_pts) / pyscan::computeTotal(m_pts) - pyscan::range_weight(rect2, b_pts) / pyscan::computeTotal(b_pts)
            << std::endl;

//        EXPECT_FLOAT_EQ(val, pyscan::range_weight(rect, m_pts) / pyscan::computeTotal(m_pts) -
//                             pyscan::range_weight(rect, b_pts) / pyscan::computeTotal(b_pts));
        EXPECT_FLOAT_EQ(rect_val, pyscan::range_weight(rect2, m_pts) / pyscan::computeTotal(m_pts) -
                                  pyscan::range_weight(rect2, b_pts) / pyscan::computeTotal(b_pts));

        EXPECT_FLOAT_EQ(val, rect_val);



        //auto [rect2, rect_val] = pyscan::max_range4<pyscan::ERectangle>(n_pts, tmp_pts, m_pts, b_pts, f);


    }

    using mx_list_t = std::vector<pyscan::MaxIntervalAlt>;

    void interval_eq(mx_list_t a1, mx_list_t a2) {
        ASSERT_EQ(a1.size(), a2.size());
        for (size_t i = 0; i < a1.size(); i++) {
            ASSERT_EQ(a1[i].get_l(), a2[i].get_l());
            ASSERT_EQ(a1[i].get_r(), a2[i].get_r());
            ASSERT_FLOAT_EQ(a1[i].get_max().get_v(), a2[i].get_max().get_v());
        }
    }

//    TEST(MaxIntervalAlt, merges) {
//        mx_list_t intervals = {
//                pyscan::MaxIntervalAlt(pyscan::Interval(7, 7, 1), pyscan::Interval(7, 7, 1), pyscan::Interval(7, 7, 1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(9, 9, -1), pyscan::Interval(9, 9, -1),
//                                       pyscan::Interval(9, 9, -1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(12, 12, 1), pyscan::Interval(12, 12, 1),
//                                       pyscan::Interval(12, 12, 1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(17, 19, 2), pyscan::Interval(17, 19, 2),
//                                       pyscan::Interval(17, 19, 2)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(36, 36, -1), pyscan::Interval(36, 36, -1),
//                                       pyscan::Interval(36, 36, -1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(38, 38, -1), pyscan::Interval(38, 38, -1),
//                                       pyscan::Interval(38, 38, -1))};
//        std::vector<size_t> s1 = {};
//        std::vector<size_t> s2 = {2, 10};
//        std::vector<size_t> s3 = {8, 10};
//        std::vector<size_t> s4 = {5, 8, 10, 39};
//        std::vector<size_t> s5 = {5, 8, 10, 37};
//
//
//        auto i1 = pyscan::reduce_merges(intervals, s1);
//        std::cout << i1 << std::endl;
//        auto i2 = pyscan::reduce_merges(intervals, s2);
//        auto i3 = pyscan::reduce_merges(intervals, s3);
//        auto i4 = pyscan::reduce_merges(intervals, s4);
//        auto i5 = pyscan::reduce_merges(intervals, s5);
//
//        mx_list_t i1_ac = {pyscan::MaxIntervalAlt(pyscan::Interval(7, 7, 1), pyscan::Interval(12, 19, 3),
//                                                  pyscan::Interval(12, 38, 1))};
//        mx_list_t i2_ac = {
//                pyscan::MaxIntervalAlt(pyscan::Interval(7, 7, 1), pyscan::Interval(7, 7, 1), pyscan::Interval(7, 9, 0)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(12, 19, 3), pyscan::Interval(12, 19, 3),
//                                       pyscan::Interval(12, 38, 1))};
//        mx_list_t i3_ac = {
//                pyscan::MaxIntervalAlt(pyscan::Interval(7, 7, 1), pyscan::Interval(7, 7, 1), pyscan::Interval(7, 7, 1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(9, 9, -1), pyscan::Interval(9, 9, -1),
//                                       pyscan::Interval(9, 9, -1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(12, 19, 3), pyscan::Interval(12, 19, 3),
//                                       pyscan::Interval(12, 38, 1))};
//
//        mx_list_t i4_ac = {
//                pyscan::MaxIntervalAlt(pyscan::Interval(7, 7, 1), pyscan::Interval(7, 7, 1), pyscan::Interval(7, 7, 1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(9, 9, -1), pyscan::Interval(9, 9, -1),
//                                       pyscan::Interval(9, 9, -1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(12, 19, 3), pyscan::Interval(12, 19, 3),
//                                       pyscan::Interval(12, 38, 1))};
//        mx_list_t i5_ac = {
//                pyscan::MaxIntervalAlt(pyscan::Interval(7, 7, 1), pyscan::Interval(7, 7, 1), pyscan::Interval(7, 7, 1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(9, 9, -1), pyscan::Interval(9, 9, -1),
//                                       pyscan::Interval(9, 9, -1)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(12, 19, 3), pyscan::Interval(12, 19, 3),
//                                       pyscan::Interval(12, 36, 2)),
//                pyscan::MaxIntervalAlt(pyscan::Interval(38, 38, -1), pyscan::Interval(38, 38, -1),
//                                       pyscan::Interval(38, 38, -1))};
////        mx_list_t  i3_ac = {pyscan::MaxIntervalAlt(pyscan::Interval(7, 7, 1), pyscan::Interval(12, 19, 3), pyscan::Interval(12, 38, 1))};
////        mx_list_t  i4_ac = {pyscan::MaxIntervalAlt(pyscan::Interval(7, 7, 1), pyscan::Interval(12, 19, 3), pyscan::Interval(12, 38, 1))};
//
//
//        interval_eq(i1, i1_ac);
//        std::cout << i2 << " || " << i2_ac << std::endl;
//        interval_eq(i2, i2_ac);
//        std::cout << i3 << " || " << i3_ac << std::endl;
//        interval_eq(i3, i3_ac);
//        std::cout << i4 << " || " << i4_ac << std::endl;
//        interval_eq(i4, i4_ac);
//        std::cout << i5 << " || " << i5_ac << std::endl;
//        interval_eq(i5, i5_ac);
//    }


    TEST(MaxIntervalAlt, merge_testing) {
        std::random_device rd;
        std::default_random_engine default_gen(rd());

        std::uniform_real_distribution<double> distr(-1.0, 1.0);

        std::vector<size_t> index;
        std::vector<double> weights;
        std::vector<pyscan::MaxIntervalAlt> curr_intervals;
        //Create a list of values
        for (size_t i = 0; i < 10; i++) {
            index.emplace_back(i);
            weights.emplace_back(distr(default_gen));
            curr_intervals.emplace_back(i, weights.back());
        }

        auto max_interval = pyscan::max_interval(index, weights);
        auto act_max_interval = pyscan::max_interval_slow(index, weights);

        double w = 0;
        for (size_t i = max_interval.get_l(); i <= max_interval.get_r(); i++) {
            w += weights[i];
        }
        EXPECT_FLOAT_EQ(w, max_interval.get_v());
        EXPECT_FLOAT_EQ(act_max_interval.get_v(), max_interval.get_v());
        EXPECT_FLOAT_EQ(act_max_interval.get_l(), max_interval.get_l());
        EXPECT_FLOAT_EQ(act_max_interval.get_r(), max_interval.get_r());

        while (curr_intervals.size() > 1) {
            std::uniform_int_distribution<size_t> range(1, curr_intervals.size() - 1);
            size_t ix = range(default_gen);
            curr_intervals[ix - 1] += curr_intervals[ix];
            curr_intervals.erase(curr_intervals.begin() + ix);
            std::cout << curr_intervals << std::endl;
        }
        double nw = 0;
        for (size_t i = curr_intervals[0].get_max().get_l(); i <= curr_intervals[0].get_max().get_r(); i++) {
            nw += weights[i];
        }
        EXPECT_FLOAT_EQ(nw, curr_intervals[0].get_max().get_v());

        EXPECT_FLOAT_EQ(max_interval.get_v(), curr_intervals[0].get_max().get_v());

        std::cout << max_interval << std::endl;
        std::cout << curr_intervals[0] << std::endl;

    }
}
