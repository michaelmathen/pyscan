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


    Point<3> cross_product(const Point<3> &p1, const Point<3> &p2) {
        return Point<3>( util::det2(p1(1), p2(1), p1(2), p2(2))
                , util::det2(p1(0), p2(0), p1(2), p2(2))
                , util::det2(p1(0), p2(0), p1(1), p2(1)), 1.0);
    }

}