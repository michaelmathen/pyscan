
#include <cmath>
#include <algorithm>
#include <vector>
#include <tuple>
#include <functional>

#include "Utilities.hpp"
#include "Statistics.hpp"
#include "HalfplaneScan.hpp"


namespace pyscan {


    template <typename T, typename Compare>
    std::vector<std::size_t> sort_permutation(const std::vector<T>& vec, Compare compare) {
        std::vector<std::size_t> p(vec.size());
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(),
                  [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
        return p;
    }

    template <typename T>
    void apply_permutation_in_place(std::vector<T>& vec, std::vector<std::size_t> const& p) {
        std::vector<bool> done(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i) {
            if (done[i])
            {
                continue;
            }
            done[i] = true;
            std::size_t prev_j = i;
            std::size_t j = p[i];
            while (i != j)
            {
                std::swap(vec[prev_j], vec[j]);
                done[j] = true;
                prev_j = j;
                j = p[j];
            }
        }
    }


    double order_function(Point<> const& p1, Point<> const& p2) {
        double y = p1[2] * p2[0] - p2[2] * p1[0];
        double x = p1[2] * p2[1] - p2[2] * p1[1];
        if (y >= 0) {
            return atan2(y, x);
        } else {
            return M_PI + atan2(y, x);
        }

    }

    std::tuple<Point<>, double> max_halfplane(
            point_list& point_net,
            point_list& red, weight_list& red_weight,
            point_list& blue, weight_list& blue_weight,
            std::function<double(double, double)> f) {
        double totalR = 0;
        double totalB = 0;
        for (auto& rw : red_weight) {
            totalR += rw;
        }
        for (auto& bw : blue_weight) {
            totalB += bw;
        }

        double max_discrepancy = -std::numeric_limits<double>::infinity();
        Point<> max_line;

        for (size_t i = 0; i < point_net.size() - 1; i++) {
            auto p_0 = point_net[i];

            auto order_f = [&p_0] (Point<> const& p1) {
                return order_function(p_0, p1);
            };
            auto pb = point_net.begin() + i + 1;

            std::sort(pb, point_net.end(), [&p_0] (Point<> const& p1, Point<> const& p2) {
                return order_f(p1) < order_f(p2);
            });

            auto l1 = intersection(p_0, *pb);

            std::vector<double> red_delta(point_net.end() - pb, 0.0);
            std::vector<double> blue_delta(point_net.end() - pb, 0.0);
            std::vector<double> angles(point_net.end() - pb, 0.0);
            for (size_t j = 0; j < angles.size(); j++) {
                angles[j] = order_f(*(pb + j));
            }
            double red_curr = 0;
            for (size_t j = 0; j < red.size(); j++) {
                auto angle_it = std::upper_bound(angles.begin(), angles.end(), order_f(red[i]));
                if (angle_it == angles.begin()) {
                    continue
                } else {
                    size_t ix = angle_it - angles.begin() - 1;
                    if (l1.above_closed(red[i])) {
                        red_delta[ix] -= red_weight[i];
                        red_curr += red_weight[i];
                    } else {
                        red_delta[ix] += red_weight[i];
                    }
                }
            }
            double blue_curr = 0;
            for (size_t j = 0; j < blue.size(); j++) {
                auto angle_it = std::upper_bound(angles.begin(), angles.end(), order_f(blue[i]));
                if (angle_it == angles.begin()) {
                    continue
                } else {
                    size_t ix = angle_it - angles.begin() - 1;
                    if (l1.above_closed(red[i])) {
                        blue_delta[ix] -= blue_weight[i];
                        blue_curr += blue_weight[i];
                    } else {
                        blue_delta[ix] += blue_weight[i];
                    }
                }
            }

            for (size_t j = 0; j < angles.size(); j++) {
                double stat = f(red_curr, blue_curr);
                if (max_discrepancy <= stat) {
                    max_line = intersection(p_0, *(j + pb));
                    max_discrepancy = stat;
                }
                red_curr += red_delta[j];
                blue_curr += blue_delta[j];
            }
        }
        return std::make_tuple(max_line, max_discrepancy);
    }


    Point<3> lift_halfspace(Point<2> const& h, Point<3> const& p) {

    }

    Point<2> dropPoint(Point<3> const& fixed_point, Point<3> const& p1) {
        /*
         * Does an affine transformation from 3 to 2.
         */
        return {p1[0] * fixed_point[2] - fixed_point[0] * p1[2],
                p1[1] * fixed_point[2] - fixed_point[1] * p1[2],
                p1[3] * fixed_point[2] - fixed_point[3] * p1[2]};
    }
}
