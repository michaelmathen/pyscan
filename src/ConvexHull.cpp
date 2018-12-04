//
// Created by mmath on 12/3/18.
//

#include "ConvexHull.hpp"


namespace pyscan {

    // TODO remove all the () so there is no division going on in here.
    inline double orientation(const pt2_t& p, const pt2_t& q, const pt2_t& r) {
        return  ((q(1) - p(1)) * (r(0) - q(0)) - (q(0) - p(0)) * (r(1) - q(1)));
//        return (q[1] * p[2] - p[1] * q[2]) * (r[0] * q[2] - q[0] * r[2]) -
//                (q[0] * p[2] - p[0] * q[2]) * (r[1] * q[2] - q[1] * r[2]);
    }


    point_list_t graham_march(point_list_t points) {
        point_list_t convex_hull;
        auto p0_iter = std::min_element(points.begin(), points.end(), [](const pt2_t &p1, const pt2_t &p2) {

            return (p1(1)  < p2(1)) || ((p1(1) == p2(1)) && (p1(0) < p2(0)));
            //return (p1[1] * p2[2] < p2[1] * p1[2]) || ((p1[1] * p2[2] == p2[1] * p1[2]) && (p1[0] * p2[2] < p2[0] * p1[2]));
        });

        auto p0 = *p0_iter;
        std::iter_swap(p0_iter, points.end() - 1);
        convex_hull.push_back(*(points.end() - 1));

        //sort by polar angle with bottom leftmost point.
        std::sort(points.begin(), points.end() - 1, [&p0](const pt2_t &p1, const pt2_t &p2) {
            double val = orientation(p0, p1, p2);
            if (util::aeq(val, 0.0)) {
                //Place the farthest point first
                return p0.square_dist(p1) > p0.square_dist(p2);
            } else {
                return val > 0;
            }
        });

        auto new_last = std::unique(points.begin(), points.end() - 1,
                                    [&p0](const pt2_t &p1, const pt2_t &p2) {
                                        //Will delete all duplicates
                                        return util::aeq(orientation(p0, p1, p2), 0.0);
                                    }
        );
        //auto new_last = points.end() - 1;
        if (new_last - points.begin() < 2) {
            std::copy(points.begin(), new_last, std::back_inserter(convex_hull));
            return convex_hull;
        }

        convex_hull.push_back(*points.begin());
        convex_hull.push_back(*(points.begin() + 1));

        for (auto it = points.begin() + 2; it != new_last; ++it) {
            while (orientation(convex_hull[convex_hull.size() - 2],
                               convex_hull[convex_hull.size() - 1],
                               *it) <= 0) {
                convex_hull.pop_back();
            }
            convex_hull.push_back(*it);
        }
        return convex_hull;
    }
}
