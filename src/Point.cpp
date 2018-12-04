//
// Created by mmath on 9/27/18.
//
#include "Point.hpp"

namespace pyscan {

    Point<2> correct_orientation(const Point<2>& pivot, const Point<2>& p) {
        double y = pivot[2] * p[0] - p[2] * pivot[0];
        if (util::alte(0.0, y)) return intersection(pivot, p);
        else return intersection(p, pivot);
    }

    Point<2> intersection(const Point<2> &p1, const Point<2> &p2) {
        return Point < 2 > {util::det2(p1[1], p1[2], p2[1], p2[2]),
                            -util::det2(p1[0], p1[2], p2[0], p2[2]),
                            util::det2(p1[0], p1[1], p2[0], p2[1])};
    }

    bool is_parallel(const Point<2> &l1, const Point<2> &l2) {
        return util::aeq(intersection(l1, l2)[2], 0.0);
    }

    std::tuple<double, double, double> normal(const Point<3> &p1, const Point<3> &p2, const Point<3> &p3) {
        auto v1 = p1 - p2;
        auto v2 = p1 - p3;

        double ax = v1(0), ay = v1(1), az = v1(2);
        double bx = v2(0), by = v2(1), bz = v2(2);
        return {util::det2(ay, az, by, bz),
                -util::det2(ax, az, bx, bz),
                util::det2(ax, ay, bx, by)};
    }

    bool crosses_segment(const Point<2> &p1, const Point<2> &p2, const Point<2> &q1, const Point<2> &q2) {
        /*
         * Checks if these two segment cross.
         */
        return intersection(q1, q2).crosses(p1, p2) && intersection(p1, p2).crosses(q1, q2);
    }

    Point<3> cross_product(const Point<3>& p1, const Point<3>& p2) {
        assert(util::aeq(p1[3], 1.0));
        assert(util::aeq(p2[3], 1.0));
        return Point <3> {util::det2(p1[1], p1[2], p2[1], p2[2]),
                            -util::det2(p1[0], p1[2], p2[0], p2[2]),
                            util::det2(p1[0], p1[1], p2[0], p2[1]),
                            1.0
        };
    }
//     double projected_distance(const Point<2>& p1, const Point<2> &p2, const Point<2> &origin, double dist) {
//        //X * line(0) + Y * line(1) + 1 = 0
//        //Y = (-1 - X * line(0)) / line(1)
//        // (-b +- sqrt(b^2 -4ac))/2a
//
//        // X^2 - 2 * origin(0) * X + origin(0) * origin(0) + Y^2 - 2*Y*origin(1) + origin(1) * origin(1) - dist*dist=0
//        //X^2 - 2 * origin(O) * X + origin(0)^2 + ((-1 - X * line(0)) / line(1))^2 - 2*(-1-X*line(0))/line(1) * origin(1) + origin(1)^2 - dist^2 = 0;
//
//        //X^2 - 2 * origin(O) * X + origin(0)^2 + (1 + X * line(0) + X^2 * line(0)^2) / line(1)^2 + 2*(1 + X*line(0))/line(1) * origin(1) + origin(1)^2 - dist^2 = 0;
//
//        //1 + X *line(0) + line(0)^2 / line(1)^2 * X^2
//
//
//        // ||p_1 - p_2||
//
//        pt2_t line = intersection(p1, p1);
//        double a = (1 + line(0) * line(0)) / (line(1) * line(1));
//        double b = (-2 * origin(0) + line(0)/(line(1) * line(1)) + 2 * origin(1) * line(0) / line(1));
//        double c = origin(0) * origin(0) + 1 / (line(0) * line(0)) + 2 * origin(1) / line(1) + origin(1) * origin(1) - dist * dist;
//
//        if (b * b - 4 * a * c <= 0) {
//            return 0;
//        }
//
//        double p1_x = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
//        double p2_x_= (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
//        // -line(0) / line(1) * px -1 = py
//        double p1_y = -line(0) / line(1) * p1_x - 1 / line(1);
//        double p2_y = -line(0) / line(1) * p2_x - 1 / line(1);
//        ret
//
//        // p1 <-> p2
//
//    }
}