//
// Created by mmath on 10/2/17.
//


//#include "../src/RectangleScan.hpp"
#include "DiskScan.hpp"
#include "Range.hpp"
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


    static auto scan = [](double m, double m_total, double b, double b_total) {
        return fabs(m / m_total - b / b_total);
    };

    TEST(DiskTest, matching) {

        const static int n_size = 50;
        const static int s_size = 1000;
        auto n_pts = pyscantest::randomPoints2(n_size);
        auto m_pts = pyscantest::randomWPoints2(s_size);
        auto b_pts = pyscantest::randomWPoints2(s_size);

        auto [d2, d2value] = pyscan::max_disk(n_pts, m_pts, b_pts, scan);
        std::cout << d2value << std::endl;
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        auto [d1, d1value] = pyscan::max_range3<pyscan::Disk>(n_pts, m_pts, b_pts, scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));

    }


    TEST(max_disk_scale, matching) {

        const static int n_size = 25;
        const static int s_size = 100;
        auto n_pts = pyscantest::randomPoints2(n_size);
        auto m_pts = pyscantest::randomWPoints2(s_size);
        auto b_pts = pyscantest::randomWPoints2(s_size);

        auto [d1, d1value] = max_disk_scale(n_pts, m_pts, b_pts, 1 / 32.0, scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
    }



//    TEST(DiskScan2, matching) {
//
//        const static int n_size = 25;
//        const static int s_size = 100;
//        auto n_pts = pyscantest::randomPoints2(n_size);
//        auto m_pts = pyscantest::randomWPoints2(s_size);
//        auto b_pts = pyscantest::randomWPoints2(s_size);
//
//        pyscan::Disk d1, d2;
//        double d1value, d2value;
//        std::tie(d1, d1value) = max_disk_lift(n_pts, m_pts, b_pts, scan);
//        std::tie(d2, d2value) = max_disk_simple(n_pts, m_pts, b_pts, scan);
//
//        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
//        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
//        EXPECT_FLOAT_EQ(d1value, d2value);
//
//    }

    void label_test_restricted(size_t n_size, size_t s_size, size_t labels) {
        auto n_lpts = pyscantest::randomLPoints2(n_size, labels);
        auto n_pts = pyscantest::removeLW(n_lpts);
        auto m_pts = pyscantest::randomLPoints2(s_size, labels);
        auto b_pts = pyscantest::randomLPoints2(s_size, labels);

        auto [d1, d1value] = pyscan::max_disk_scale_labeled(n_lpts, m_pts, b_pts, true, 0.2, scan);
        auto [d2, d2value] = pyscan::max_disk_scale_labeled_alt(n_pts, m_pts, b_pts, 0.2, scan);

        auto [d3, d3value] = pyscan::max_disk_scale_slow_labeled(n_pts, m_pts, b_pts, 0.2, 0.4, scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d3value, evaluate_range(d3, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d2value, d1value);
        EXPECT_FLOAT_EQ(d3value, d1value);
        EXPECT_FLOAT_EQ(d3value, d2value);
    }


    TEST(max_disk_scale_labeled, matching) {
       // label_test_restricted(100, 10000, 100);

    }

    void label_test(size_t n_size, size_t s_size, size_t labels) {
      auto n_pts = pyscantest::randomPoints2(n_size);
      auto m_pts = pyscantest::randomLPoints2(s_size, labels);
      auto b_pts = pyscantest::randomLPoints2(s_size, labels);
      //std::cout << m_pts << std::endl;
      auto [d1, d1value] = pyscan::max_disk_labeled(n_pts, m_pts, b_pts, scan);
      EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
      auto [d2, d2value] = pyscan::max_range3<pyscan::Disk, pyscan::LPoint>(n_pts, m_pts, b_pts, scan);

      std::cout << d2value << " " << d1value << std::endl;
      EXPECT_FLOAT_EQ(d1value, d2value);

      EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
    }


//    void label_test2(size_t n_size, size_t s_size, size_t labels) {
//        auto n_pts = pyscantest::randomPoints2(n_size);
//        auto m_pts = pyscantest::randomLPoints2(s_size, labels);
//        auto b_pts = pyscantest::randomLPoints2(s_size, labels);
//
//        pyscan::Disk d1, d2;
//        double d1value, d2value;
//        std::tie(d1, d1value) = pyscan::max_disk_lift_labeled(n_pts, m_pts, b_pts, scan);
//        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
//
//        std::tie(d2, d2value) = pyscan::max_disk_simple_labeled(n_pts, m_pts, b_pts, scan);
//        EXPECT_FLOAT_EQ(d1value, d2value);
//
//        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
//    }


    TEST(DiskTestLabel, special) {
        pyscan::point_list_t n_pts = {pyscan::Point<2>(0.55158, 0.553599, 1.0), pyscan::Point<2>(0.551583, 0.553599, 1.0), pyscan::Point<2>(0.823782, 0.630888, 1.0), pyscan::Point<2>(0.823947, 0.631124, 1.0), pyscan::Point<2>(0.824017, 0.631223, 1.0),
                                        pyscan::Point<2>(0.824113, 0.631359, 1.0), pyscan::Point<2>(0.824253, 0.631558, 1.0), pyscan::Point<2>(0.824279, 0.631595, 1.0), pyscan::Point<2>(0.824445, 0.631831, 1.0), pyscan::Point<2>(0.824489, 0.631893, 1.0),
                                        pyscan::Point<2>(0.82461, 0.632066, 1.0), pyscan::Point<2>(0.824724, 0.632228, 1.0), pyscan::Point<2>(0.82476, 0.63228, 1.0), pyscan::Point<2>(0.824765, 0.632302, 1.0), pyscan::Point<2>(0.824818, 0.632538, 1.0),
                                        pyscan::Point<2>(0.824871, 0.632773, 1.0), pyscan::Point<2>(0.824924, 0.633009, 1.0), pyscan::Point<2>(0.82496, 0.633169, 1.0), pyscan::Point<2>(0.824977, 0.633245, 1.0), pyscan::Point<2>(0.824994, 0.633322, 1.0),
                                        pyscan::Point<2>(0.477515, 0.914504, 1.0), pyscan::Point<2>(0.477521, 0.914494, 1.0), pyscan::Point<2>(0.682031, 0.964714, 1.0), pyscan::Point<2>(0.682032, 0.964716, 1.0),
                                        pyscan::Point<2>(0.411749, 0.62699, 1.0), pyscan::Point<2>(0.411759, 0.626994, 1.0), pyscan::Point<2>(0.374895, 0.611139, 1.0), pyscan::Point<2>(0.374899, 0.611137, 1.0), pyscan::Point<2>(0.688543, 0.582471, 1.0),
                                        pyscan::Point<2>(0.688556, 0.582447, 1.0), pyscan::Point<2>(0.682023, 0.964699, 1.0),pyscan::Point<2>(0.682028, 0.964697, 1.0), pyscan::Point<2>(0.595725, 0.63875, 1.0), pyscan::Point<2>(0.595727, 0.638749, 1.0),
                                        pyscan::Point<2>(0.694211, 0.965976, 1.0), pyscan::Point<2>(0.694212, 0.965976, 1.0), pyscan::Point<2>(0.538954, 0.35394, 1.0), pyscan::Point<2>(0.538961, 0.353934, 1.0), pyscan::Point<2>(0.348578, 0.808888, 1.0),
                                        pyscan::Point<2>(0.348579, 0.808889, 1.0), pyscan::Point<2>(0.63142, 0.423662, 1.0), pyscan::Point<2>(0.631423, 0.423664, 1.0), pyscan::Point<2>(0.553365, 0.920366, 1.0), pyscan::Point<2>(0.553366, 0.920366, 1.0),
                                        pyscan::Point<2>(0.717556, 0.25, 1.0), pyscan::Point<2>(0.717567, 0.24999, 1.0), pyscan::Point<2>(0.717651,0.249754, 1.0), pyscan::Point<2>(0.717792, 0.249784, 1.0), pyscan::Point<2>(0.717798, 0.249779, 1.0),
                                        pyscan::Point<2>(0.717967, 0.249937, 1.0), pyscan::Point<2>(0.358201, 0.632762, 1.0), pyscan::Point<2>(0.358201, 0.632762, 1.0), pyscan::Point<2>(0.54108, 0.350768, 1.0), pyscan::Point<2>(0.541133, 0.351003, 1.0),
                                        pyscan::Point<2>(0.541187, 0.351239, 1.0), pyscan::Point<2>(0.541195, 0.351272, 1.0), pyscan::Point<2>(0.541315, 0.351265, 1.0), pyscan::Point<2>(0.541551, 0.35125, 1.0), pyscan::Point<2>(0.541729, 0.351239, 1.0),
                                        pyscan::Point<2>(0.541787, 0.351236, 1.0), pyscan::Point<2>(0.542022, 0.351221, 1.0), pyscan::Point<2>(0.542162, 0.351212, 1.0),
                                        pyscan::Point<2>(0.646433, 0.490692, 1.0), pyscan::Point<2>(0.646669, 0.490773, 1.0), pyscan::Point<2>(0.64687, 0.490842, 1.0), pyscan::Point<2>(0.646905, 0.490768, 1.0), pyscan::Point<2>(0.646972, 0.49062, 1.0),
                                        pyscan::Point<2>(0.64708, 0.490384, 1.0), pyscan::Point<2>(0.64714, 0.490254, 1.0), pyscan::Point<2>(0.647188, 0.490149, 1.0), pyscan::Point<2>(0.647376, 0.490285, 1.0), pyscan::Point<2>(0.647513, 0.490384, 1.0),
                                        pyscan::Point<2>(0.647612, 0.490456, 1.0), pyscan::Point<2>(0.647838, 0.49062, 1.0), pyscan::Point<2>(0.647847, 0.490627, 1.0), pyscan::Point<2>(0.647936, 0.490692, 1.0), pyscan::Point<2>(0.534208, 0.346251, 1.0),
                                        pyscan::Point<2>(0.534214, 0.346257, 1.0), pyscan::Point<2>(0.55882, 0.6476, 1.0), pyscan::Point<2>(0.558824, 0.647612, 1.0), pyscan::Point<2>(0.34864, 0.808983, 1.0), pyscan::Point<2>(0.348641, 0.808986, 1.0),
                                        pyscan::Point<2>(0.373765, 0.610296, 1.0), pyscan::Point<2>(0.373767, 0.610296, 1.0), pyscan::Point<2>(0.211007, 0.177877, 1.0), pyscan::Point<2>(0.211007, 0.177876, 1.0), pyscan::Point<2>(0.597078, 0.636954, 1.0), pyscan::Point<2>(0.597079, 0.636954, 1.0),
                                        pyscan::Point<2>(0.414363, 0.632891, 1.0), pyscan::Point<2>(0.414364, 0.632896, 1.0), pyscan::Point<2>(0.633982, 0.66672, 1.0), pyscan::Point<2>(0.634218, 0.666723, 1.0), pyscan::Point<2>(0.634229, 0.666723, 1.0),
                                        pyscan::Point<2>(0.634353, 0.666488, 1.0), pyscan::Point<2>(0.634453, 0.666298, 1.0), pyscan::Point<2>(0.634478, 0.666253, 1.0), pyscan::Point<2>(0.634602, 0.666017, 1.0), pyscan::Point<2>(0.634689, 0.665853, 1.0),
                                        pyscan::Point<2>(0.634737, 0.665781, 1.0), pyscan::Point<2>(0.634925, 0.665616, 1.0), pyscan::Point<2>(0.635004, 0.665545, 1.0), pyscan::Point<2>(0.635161, 0.665408, 1.0), pyscan::Point<2>(0.635272, 0.66531, 1.0),
                                        pyscan::Point<2>(0.635396, 0.665201, 1.0), pyscan::Point<2>(0.63554, 0.665074, 1.0), pyscan::Point<2>(0.635632, 0.664993, 1.0), pyscan::Point<2>(0.635808, 0.664838, 1.0), pyscan::Point<2>(0.635868, 0.664786, 1.0),
                                        pyscan::Point<2>(0.635979, 0.664688, 1.0), pyscan::Point<2>(0.636103, 0.664675, 1.0), pyscan::Point<2>(0.636339, 0.66465, 1.0), pyscan::Point<2>(0.636575, 0.664625, 1.0), pyscan::Point<2>(0.636793, 0.664603, 1.0),
                                        pyscan::Point<2>(0.63681, 0.664615, 1.0), pyscan::Point<2>(0.637046, 0.664779, 1.0), pyscan::Point<2>(0.637132, 0.664838, 1.0), pyscan::Point<2>(0.637282, 0.664943, 1.0), pyscan::Point<2>(0.637471, 0.665074, 1.0), pyscan::Point<2>(0.637518, 0.665107, 1.0),
                                        pyscan::Point<2>(0.637753, 0.665271, 1.0), pyscan::Point<2>(0.637809, 0.66531, 1.0), pyscan::Point<2>(0.637989, 0.665435, 1.0), pyscan::Point<2>(0.638073, 0.665493, 1.0), pyscan::Point<2>(0.638225, 0.665417, 1.0),
                                        pyscan::Point<2>(0.638438, 0.66531, 1.0), pyscan::Point<2>(0.63846, 0.665299, 1.0), pyscan::Point<2>(0.638696, 0.66518, 1.0), pyscan::Point<2>(0.638779, 0.665138, 1.0), pyscan::Point<2>(0.677457, 0.491775, 1.0),
                                        pyscan::Point<2>(0.677472, 0.491767, 1.0), pyscan::Point<2>(0.494189, 0.639286, 1.0), pyscan::Point<2>(0.494288, 0.639492, 1.0), pyscan::Point<2>(0.494425, 0.639474, 1.0), pyscan::Point<2>(0.494661, 0.639442, 1.0),
                                        pyscan::Point<2>(0.494896, 0.63941, 1.0), pyscan::Point<2>(0.494954, 0.640344, 1.0), pyscan::Point<2>(0.495046, 0.640229, 1.0), pyscan::Point<2>(0.495114, 0.639892, 1.0), pyscan::Point<2>(0.495132, 0.639378, 1.0),
                                        pyscan::Point<2>(0.495132, 0.63983, 1.0), pyscan::Point<2>(0.495132, 0.639378, 1.0), pyscan::Point<2>(0.495132, 0.640121, 1.0), pyscan::Point<2>(0.495132, 0.640371, 1.0), pyscan::Point<2>(0.495132, 0.640121, 1.0), pyscan::Point<2>(0.495132, 0.640371, 1.0),
                                        pyscan::Point<2>(0.495132, 0.640121, 1.0), pyscan::Point<2>(0.495132, 0.640371, 1.0), pyscan::Point<2>(0.495132, 0.640121, 1.0), pyscan::Point<2>(0.495152, 0.639758, 1.0), pyscan::Point<2>(0.495218, 0.639522, 1.0),
                                        pyscan::Point<2>(0.495234, 0.639993, 1.0), pyscan::Point<2>(0.495264, 0.63936, 1.0), pyscan::Point<2>(0.495304, 0.639906, 1.0), pyscan::Point<2>(0.495368, 0.640406, 1.0), pyscan::Point<2>(0.495372, 0.640407, 1.0)};


        auto m_pts = pyscantest::randomLPoints2(1000, 50);
        auto b_pts = pyscantest::randomLPoints2(1000, 50);
        //std::cout << m_pts << std::endl;
        auto [d1, d1value] = pyscan::max_disk_labeled(n_pts, m_pts, b_pts, scan);
        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        auto [d2, d2value] = pyscan::max_range3<pyscan::Disk, pyscan::LPoint>(n_pts, m_pts, b_pts, scan);

        std::cout << d2value << " " << d1value << std::endl;
        EXPECT_FLOAT_EQ(d1value, d2value);

        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
    }

    TEST(DiskTestLabel, allsamelabel) {
      label_test(25, 100, 1);
    }

    TEST(DiskTestLabel, mixedlabels) {
      label_test(25, 100, 5);
      label_test(25, 1000, 100);
    }


//    TEST(max_rdisk_lifted, mixedlabels) {
//        label_test_restrictedt(25, 100, 10);
//    }


//    TEST(DiskScan2, allsamelabel) {
//        label_test2(25, 100, 1);
//    }
//
//    TEST(DiskScan2, mixedlabels) {
//        label_test2(25, 100, 5);
//    }

    TEST(DiskTestLabel, alluniquelabel) {

        const static int n_size = 25;
        const static int s_size = 100;

        auto n_pts = pyscantest::randomPoints2(n_size);
        auto m_lpts = pyscantest::randomLPointsUnique2(s_size);
        auto b_lpts = pyscantest::randomLPointsUnique2(s_size);

        auto m_pts = pyscantest::removeLabels(m_lpts);
        auto b_pts = pyscantest::removeLabels(b_lpts);

        auto [d1, d1value] = pyscan::max_disk_labeled(n_pts, m_lpts, b_lpts, scan);
        auto [d2, d2value] = pyscan::max_disk(n_pts, m_pts, b_pts, scan);

        auto [d3, d3value] = pyscan::max_range3<pyscan::Disk, pyscan::LPoint>(n_pts, m_lpts, b_lpts, scan);
        auto [d4, d4value] = pyscan::max_range3<pyscan::Disk, pyscan::WPoint>(n_pts, m_pts, b_pts, scan);


        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d3value, evaluate_range(d3, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d4value, evaluate_range(d4, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d1value, d2value);
        EXPECT_FLOAT_EQ(d1value, d3value);
        EXPECT_FLOAT_EQ(d2value, d4value);
        EXPECT_FLOAT_EQ(d3value, d4value);

    }

    TEST(DiskTestLabelSlow, matching) {

        const static int n_size = 25;
        const static int s_size = 100;
        auto n_pts = pyscantest::randomPoints2(n_size);
        auto m_lpts = pyscantest::randomLPointsUnique2(s_size);
        auto b_lpts = pyscantest::randomLPointsUnique2(s_size);

        auto m_pts = pyscantest::removeLabels(m_lpts);
        auto b_pts = pyscantest::removeLabels(b_lpts);

        pyscan::Disk d1, d2;
        double d1value, d2value;
        std::tie(d1, d1value) = pyscan::max_range3<pyscan::Disk, pyscan::LPoint>(n_pts, m_lpts, b_lpts, scan);
        std::tie(d2, d2value) = pyscan::max_range3<pyscan::Disk, pyscan::WPoint>(n_pts, m_pts, b_pts, scan);

        EXPECT_FLOAT_EQ(d1value, evaluate_range(d1, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d2value, evaluate_range(d2, m_pts, b_pts, scan));
        EXPECT_FLOAT_EQ(d1value, d2value);

    }


}
