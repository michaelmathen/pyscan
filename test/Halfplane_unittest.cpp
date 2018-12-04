
#include "Test_Utilities.hpp"

#include <limits.h>
#include <random>
#include <iostream>

#include "Range.hpp"
#include "HalfSpaceScan.hpp"

#include "gtest/gtest.h"
namespace {


    TEST(halfspace3_t, lifting) {
        auto net = pyscantest::randomPoints3(100);
        auto pivot = pyscantest::randomPoints3(1);
        pyscan::point_list_t pts;
        for (auto p : net) {
            pts.push_back(pyscan::drop_point(p, pivot[0]));
        }

        for (size_t i = 0; i < net.size() - 1; i++) {
            for (size_t j = i + 1; j < net.size(); j++) {
                pyscan::halfspace2_t h(pts[i], pts[j]);
                pyscan::halfspace3_t h_lift = pyscan::halfspace3_t(pyscan::lift_half_space(h, pivot[0]));

                EXPECT_NEAR(h_lift.get_coords().evaluate(net[i]), 0.0, 10e-13);
                EXPECT_NEAR(h_lift.get_coords().evaluate(net[j]), 0.0, 10e-13);
                EXPECT_NEAR(h_lift.get_coords().evaluate(pivot[0]), 0.0, 10e-13);
                //std::cout << h << std::endl;

                for (size_t k = 0; k < net.size(); k++) {
                    //std::cout << pts[k] << net[k] << std::endl;
                    //std::cout << h.contains(pts[k].flip_orientation()) << " " << h_lift.contains(net[k]) << std::endl;
                    if (k != i && k != j) {
                        //std::cout << pivot[0](2) - net[k](2) << std::endl;
                        EXPECT_EQ(h.contains(pts[k].flip_orientation()), h_lift.contains(net[k]));
                   }
                }
            }
        }
    }
    pyscan::discrepancy_func_t stat = [](double m, double m_total, double b, double b_total) {
        return std::abs(m / m_total - b / b_total);
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

        auto [fast_h, fast_value] = pyscan::max_halfspace(n_pts, m_pts, b_pts, stat);
        auto [slow_h, slow_value] = pyscan::max_halfspace_simple(n_pts, m_pts, b_pts, stat);

        EXPECT_FLOAT_EQ(fast_value, pyscan::evaluate_range(fast_h, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(slow_value, pyscan::evaluate_range(slow_h, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(fast_value, slow_value);
    }


      TEST(max_halfplane_labeled, discrepancy) {

        const static int n_size = 25;
        const static int s_size = 100;
        for (int i = 0; i < 40; i++) {
            auto n_pts = pyscantest::randomPoints2(n_size);
            auto m_pts = pyscantest::randomLPoints2(s_size, 10);
            auto b_pts = pyscantest::randomLPoints2(s_size, 10);



            auto [fast_h, fast_value] = pyscan::max_halfplane_labeled(n_pts, m_pts, b_pts, stat);
            auto [slow_h, slow_value] = pyscan::max_halfplane_labeled_simple(n_pts, m_pts, b_pts, stat);

            EXPECT_FLOAT_EQ(slow_value, pyscan::evaluate_range(slow_h, m_pts, b_pts, stat));
            EXPECT_FLOAT_EQ(fast_value, pyscan::evaluate_range(fast_h, m_pts, b_pts, stat));
            EXPECT_FLOAT_EQ(fast_value, slow_value);
        }

    }

    void test_halfspace_labels(size_t n_size, size_t s_size, size_t num_labels) {

        auto n_pts = pyscantest::randomPoints3(n_size);
        auto m_pts = pyscantest::randomLPoints3(s_size, num_labels);
        auto b_pts = pyscantest::randomLPoints3(s_size, num_labels);


        auto [fast_h, fast_value] = pyscan::max_halfspace_labeled(n_pts, m_pts, b_pts, stat);
        auto [slow_h, slow_value] = pyscan::max_halfspace_labeled_simple(n_pts, m_pts, b_pts, stat);

        std::cout << fast_h.get_coords() << " " << fast_value << std::endl;
        std::cout << slow_h.get_coords() << " " << slow_value << std::endl;

        EXPECT_FLOAT_EQ(slow_value, pyscan::evaluate_range(slow_h, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(fast_value, pyscan::evaluate_range(fast_h, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(fast_value, slow_value);

    }

    TEST(max_halfspace_labeled, mixed_label) {
        test_halfspace_labels(25, 100, 10);

    }

    TEST(max_halfspace_labeled, same_label) {
        test_halfspace_labels(25, 100, 1);
    }

    TEST(max_halfspace_labeled, unique_label) {
        size_t n_size = 25;
        size_t s_size = 100;

        auto n_pts = pyscantest::randomPoints3(n_size);
        auto m_pts = pyscantest::randomLPointsUnique3(s_size);
        auto b_pts = pyscantest::randomLPointsUnique3(s_size);


        auto [fast_h, fast_value] = pyscan::max_halfspace_labeled(n_pts, m_pts, b_pts, stat);
        auto [slow_h, slow_value] = pyscan::max_halfspace_labeled_simple(n_pts, m_pts, b_pts, stat);

        std::cout << fast_h.get_coords() << " " << fast_value << std::endl;
        std::cout << slow_h.get_coords() << " " << slow_value << std::endl;

        EXPECT_FLOAT_EQ(slow_value, pyscan::evaluate_range(slow_h, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(fast_value, pyscan::evaluate_range(fast_h, m_pts, b_pts, stat));
        EXPECT_FLOAT_EQ(fast_value, slow_value);
    }


    inline pyscan::pt3_t lift_pt(const pyscan::pt2_t &pt) {
        double x = pt(0), y = pt(1);
        return pyscan::pt3_t(x, y, x * x + y * y, 1.0);
    }

    TEST(Remap, remap) {
        auto n_pts = pyscantest::randomPoints2(4);

        pyscan::pt3_t lifted_pt1 = lift_pt(n_pts[0]);
        pyscan::pt3_t lifted_pt2 = lift_pt(n_pts[1]);
        pyscan::pt3_t lifted_pt3 = lift_pt(n_pts[2]);
        pyscan::pt3_t lifted_pt4 = lift_pt(n_pts[3]);

        auto new_x = pyscan::pt3_t (-2 * lifted_pt1(0), -2 * lifted_pt1(1), 1.0, 1.0).normalize();
        auto new_z = cross_product(new_x, pyscan::pt3_t(0.0, 1.0, 0.0, 1.0)).normalize();
        auto new_y = cross_product(new_x, new_z).normalize();

        auto proj = [&] (pyscan::pt3_t const& pt) {
            return pyscan::pt3_t(new_x.pdot(pt), new_y.pdot(pt), new_z.pdot(pt), 1.0);
        };

        //auto s_pts = pyscantest::randomWPoints3(100);
        //auto s_pts = pyscantest::randomWPoints3(100);

        pyscan::halfspace3_t proj_h(proj(lifted_pt2), proj(lifted_pt3), proj(lifted_pt4));

        pyscan::halfspace3_t h(pyscan::Point<3>(
                proj_h[0] * new_x[0] + proj_h[1] * new_y[0] + proj_h[2] * new_z[0],
                proj_h[0] * new_x[1] + proj_h[1] * new_y[1] + proj_h[2] * new_z[1],
                proj_h[0] * new_x[2] + proj_h[1] * new_y[2] + proj_h[2] * new_z[2],
                proj_h[3]));

        pyscan::halfspace3_t alt_h(lifted_pt2, lifted_pt3, lifted_pt4);

        std::cout << alt_h.get_coords() << std::endl;
        std::cout << h.get_coords() << std::endl;
    }


//    TEST(halfSpace3_t, orientation) {
//
//        auto net = pyscantest::randomPoints3(3);
//        auto m_pts = pyscantest::randomWPoints3(100);
//        auto b_pts = pyscantest::randomWPoints3(100);
//        pyscan::halfspace3_t plane(net[0], net[1], net[2]);
//        auto upside_plane = pyscan::halfspace3_t(plane.get_coords().orient_up(2));
//        std::cout << upside_plane.get_coords() << " " << plane.get_coords() << std::endl;
//        EXPECT_FLOAT_EQ(pyscan::evaluate_range(plane, m_pts, b_pts, stat), pyscan::evaluate_range(upside_plane, m_pts, b_pts, stat));
//    }
//
//    TEST(halfSpace3_t, evaluate_range_labels) {
//
//        auto net = pyscantest::randomPoints3(3);
//        auto m_pts = pyscantest::randomLPoints3(100, 10);
//        auto b_pts = pyscantest::randomLPoints3(100, 10);
//        pyscan::halfspace3_t plane(net[0], net[1], net[2]);
//        auto upside_plane = pyscan::halfspace3_t(plane.get_coords().orient_up(2));
//        std::cout << upside_plane.get_coords() << " " << plane.get_coords() << std::endl;
//        EXPECT_FLOAT_EQ(pyscan::evaluate_range(plane, m_pts, b_pts, stat), pyscan::evaluate_range(upside_plane, m_pts, b_pts, stat));
//    }

}