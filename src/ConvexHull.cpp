//
// Created by mmath on 12/3/18.
//

#include "ConvexHull.hpp"


namespace pyscan {

    // TODO remove all the () so there is no division going on in here.
    inline double orientation(const pt2_t& p, const pt2_t& q, const pt2_t& r) {
        return  ((q(1) - p(1)) * (r(0) - q(0)) - (q(0) - p(0)) * (r(1) - q(1)));
    }


    point_list_t graham_march(point_list_t points) {
        if (points.size() < 2) {
            return points;
        }

        std::sort(points.begin(), points.end(), [](const pt2_t &p1, const pt2_t &p2) {
            return (p1(0) < p2(0)) || ((p1(0) == p2(0)) && (p1(1) < p2(1)));
        });

        //auto new_last = points.end() - 1;

        auto keep_left = [](point_list_t &points) {
            point_list_t convex_hull;
            for (auto& p : points) {
                while (convex_hull.size() > 1 &&
                       orientation(convex_hull[convex_hull.size() - 2],
                                   convex_hull[convex_hull.size() - 1],
                                   p) <= 0) {
                    convex_hull.pop_back();
                }
                convex_hull.push_back(p);
            }
            return convex_hull;
        };
        auto upper_hull = keep_left(points);
        std::reverse(points.begin(), points.end());
        auto lower_hull = keep_left(points);
        std::copy(lower_hull.begin() + 1, lower_hull.end() - 1, std::back_inserter(upper_hull));
        return upper_hull;
    }
}
