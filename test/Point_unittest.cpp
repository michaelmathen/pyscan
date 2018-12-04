//
// Created by mmath on 5/22/18.
//

#include "Point.hpp"
#include "Segment.hpp"
#include "Test_Utilities.hpp"

#include <limits.h>
#include <random>
#include <iostream>



#include "gtest/gtest.h"

namespace {


    TEST(Point, intersection) {
        pyscan::Point<> p1(1.0, 1.0, 1.0);
        pyscan::Point<> p2(1.0, 0.0, 1.0);
        pyscan::Point<> h1(1.0, 0.0, -1.0);
        EXPECT_TRUE(pyscan::intersection(p1, p2).approx_eq(h1));

        pyscan::Point<> p3(0.0, 1.0, 0.0);
        pyscan::Point<> p4(1.0, 0.0, 1.0);
        pyscan::Point<> h2(1.0, 0.0, -1.0);
        EXPECT_TRUE(pyscan::intersection(p3, p4).approx_eq(h2));

        pyscan::Point<> p5(2.0, 1.0, -1.0);
        pyscan::Point<> p6(1.0, 1.0, 3.0);
        pyscan::Point<> h3(4.0, -7.0, 1.0);
        EXPECT_TRUE(pyscan::intersection(p5, p6).approx_eq(h3));
    }

    TEST(Point, dot) {
        pyscan::Point<> p1(1.0, -1.0, 3.0);
        pyscan::Point<> p2(2.0, 0.0, 1.0);
        EXPECT_FLOAT_EQ(p1.evaluate(p2), 5.0);
        pyscan::Point<> p3(1.0, -1.0, 3.0);
        pyscan::Point<> p4(0.0, 2.0, 0.0);
        EXPECT_FLOAT_EQ(p3.evaluate(p4), -2.0);
    }

    TEST(Point, parallel) {
        pyscan::Point<> p1(0.0, 1.0, 1.0);
        pyscan::Point<> p2(0.0, 1.0, 0.0);
        EXPECT_TRUE(pyscan::is_parallel(p1, p2));
        pyscan::Point<> p5(0.0, 1.0, 0.0);
        pyscan::Point<> p6(0.0, 1.0, 0.0);
        EXPECT_TRUE(pyscan::is_parallel(p5, p6));
        pyscan::Point<> p3(0.0, 1.0, 1.0);
        pyscan::Point<> p4(-1.0, 1.0, 1.0);
        EXPECT_TRUE(!pyscan::is_parallel(p3, p4));
    }

    TEST(Point, above) {
        // Point is oriented down
        pyscan::Point<> h1(-1.0, -1.0, 0.0);
        pyscan::Point<> p1(0.0, -1.0, 2.0);
        pyscan::Point<> p2(0.0, 0.0, 1.0);
        pyscan::Point<> p3(-1.0, 1.0, 1.0);
        pyscan::Point<> p4(-1.0, 2.0, 1.0);
        EXPECT_TRUE(h1.above(p1));
        EXPECT_TRUE(!h1.above(p2));
        EXPECT_TRUE(!h1.above(p3));
        EXPECT_TRUE(!h1.above(p4));
    }

    TEST(Point, above_flipped) {
        pyscan::Point<> h1(-1.0, -1.0, 0.0);
        pyscan::Point<> p1(0.0, 1.0, -2.0);
        pyscan::Point<> p2(0.0, 0.0, -1.0);
        pyscan::Point<> p3(1.0, -1.0, -1.0);
        pyscan::Point<> p4(1.0, -2.0, -1.0);
        EXPECT_TRUE(!h1.above(p1));
        EXPECT_TRUE(!h1.above(p2));
        EXPECT_TRUE(!h1.above(p3));
        EXPECT_TRUE(h1.above(p4));
    }

    TEST(Point, above_closed) {
        pyscan::Point<> h1(-1.0, -1.0, 0.0);
        pyscan::Point<> p1(0.0, -1.0, 2.0);
        pyscan::Point<> p2(0.0, 0.0, 1.0);
        pyscan::Point<> p3(-1.0, 1.0, 1.0);
        pyscan::Point<> p4(-1.0, 2.0, 1.0);
        EXPECT_TRUE(h1.above_closed(p1));
        EXPECT_TRUE(h1.above_closed(p2));
        EXPECT_TRUE(h1.above_closed(p3));
        EXPECT_TRUE(!h1.above_closed(p4));

        //Infinite point on line
        pyscan::Point<> p5(-1.0, 1.0, 0.0);
        //Infinite point on line.
        pyscan::Point<> p6(1.0, -1.0, 0.0);
        //Infinite point below line
        pyscan::Point<> p7(1.0, -2.0, 0.0);
        pyscan::Point<> p8(-2.0, 1.0, 0.0);

        //Points above the line
        pyscan::Point<> p9(-1.0, 2.0, 0.0);
        pyscan::Point<> p10(1.0, 1.0, 0.0);

        EXPECT_TRUE(h1.above_closed(p5));
        EXPECT_TRUE(h1.above_closed(p6));

        EXPECT_TRUE(h1.above_closed(p7));
        EXPECT_TRUE(h1.above_closed(p8));
        EXPECT_TRUE(!h1.above_closed(p9));
        EXPECT_TRUE(!h1.above_closed(p10));
    }


    TEST(Point, below_closed) {
        pyscan::Point<> h1(-1.0, -1.0, 0.0);
        pyscan::Point<> p1(0.0, -1.0, 2.0);
        pyscan::Point<> p2(0.0, 0.0, 1.0);
        pyscan::Point<> p3(-1.0, 1.0, 1.0);
        pyscan::Point<> p4(-1.0, 2.0, 1.0);
        EXPECT_TRUE(!h1.below_closed(p1));
        EXPECT_TRUE(h1.below_closed(p2));
        EXPECT_TRUE(h1.below_closed(p3));
        EXPECT_TRUE(h1.below_closed(p4));
    }

    TEST(Point, parallel_lte) {
        //All of these lines are parallel to each other, but have different normals.
        pyscan::Point<> h1(1.0, 1.0, 1.0);
        pyscan::Point<> h2(1.0, 1.0, 0.0);
        pyscan::Point<> h3(2.0, 2.0, 0.0);
        pyscan::Point<> h4(-2.0, -2.0, -2.0);

        pyscan::Point<> h5(2.0, 2.0, 2.0);

        //Same line
        EXPECT_TRUE(!h2.parallel_lt(h3));

        // h1 passes below h2
        EXPECT_TRUE(h1.parallel_lt(h2));
        EXPECT_TRUE(!h2.parallel_lt(h1));

        //h3 is the same line as h2, but just a different scale
        EXPECT_TRUE(h1.parallel_lt(h3));
        EXPECT_TRUE(!h3.parallel_lt(h1));

        //h5 is the same line as h1
        EXPECT_TRUE(h5.parallel_lt(h2));
        EXPECT_TRUE(!h2.parallel_lt(h5));
        EXPECT_TRUE(h5.parallel_lt(h3));
        EXPECT_TRUE(!h3.parallel_lt(h5));

        //h4 has an opposite direction as h1, but is the same line
        EXPECT_TRUE(!h4.parallel_lt(h2));
        EXPECT_TRUE(!h2.parallel_lt(h4));
    }

    TEST(POINT, orientation) {
        auto net = pyscantest::randomPoints3(100);
        auto pivot = pyscantest::randomPoints3(1);
        for (auto p : net) {
            EXPECT_FLOAT_EQ(pivot[0].evaluate(p), -pivot[0].evaluate(p.flip_orientation()));
        }
    }
//    TEST(Point, above_closed_interval) {
//        pyscan::Point<> h1(1.0, 1.0, 0.0);
//        pyscan::Point<> h2(1.0, 1.0, -1.0);
//        pyscan::Point<> p1(0.0, -1.0, 1.0);
//        pyscan::Point<> p2(2.0, 1.0, 1.0);
//        EXPECT_TRUE(!pyscan::above_closed_interval(h1, h2, p1, p2));
//        EXPECT_TRUE(pyscan::above_closed_interval(h2, h1, p1, p2));
//
//        pyscan::Point<> h3(1.0, 2.0, 0.0);
//        pyscan::Point<> h4(1.0, 1.0, 0.0);
//        pyscan::Point<> p3(0.0, 0.0, 1.0);
//        pyscan::Point<> p4(1.0, -1.0, 0.0);
//        pyscan::Point<> p5(-2.0, 1.0, 0.0);
//        EXPECT_TRUE(pyscan::above_closed_interval(h3, h4, p3, p4));
//        EXPECT_TRUE(pyscan::above_closed_interval(h4, h3, p4, p5));
//    }

}