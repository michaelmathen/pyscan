
#include <cmath>
#include <algorithm>
#include <vector>
#include <tuple>
#include <functional>
#include <cassert>
#include <iostream>

#include "HalfplaneScan.hpp"

//#define DEBUG 

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
            while (i != j) {
                std::swap(vec[prev_j], vec[j]);
                done[j] = true;
                prev_j = j;
                j = p[j];
            }
        }
    }



    std::tuple<Point<>, double> max_halfplane(
            point_list& point_net,
            point_list& red,
            weight_list& red_w,
            point_list& blue,
            weight_list& blue_w,
            std::function<double(double, double)> const& f) {

        /*
         *
         */
        assert(red.size() == red_w.size());
        assert(blue.size() == blue_w.size());

        double red_total = std::accumulate(red_w.begin(), red_w.end(), 0.0, std::plus<>());
        double blue_total = std::accumulate(blue_w.begin(), blue_w.end(), 0.0, std::plus<>());;

        double max_discrepancy = -std::numeric_limits<double>::infinity();
        Pt2 max_line;

        for (auto p_it = point_net.begin(); p_it != point_net.end() - 1; p_it++) {
            auto p_0 = *p_it;

            auto order_f = [&] (Point<> const& p1) {
                return order_function(p_0, p1);
            };

            auto pb = p_it + 1;
            std::sort(pb, point_net.end(), [&] (Point<> const& p1, Point<> const& p2) {
                return order_f(p1) < order_f(p2);
            });

            auto l1 = correct_orientation(p_0, *pb);
            std::vector<double> red_delta(point_net.end() - pb, 0.0);
            std::vector<double> blue_delta(point_net.end() - pb, 0.0);
            std::vector<double> angles(point_net.end() - pb, 0.0);
            for (size_t j = 0; j < angles.size(); j++) {
                angles[j] = order_f(*(pb + j));
            }
            auto set_delta = [&] (auto& deltas, auto b_it, auto e_it, auto w_it) {
                double w_curr = 0;
                for (; b_it != e_it; ++b_it) {
                    auto angle_it = std::upper_bound(angles.begin(), angles.end(), order_f(*b_it));
                    if (angle_it == angles.begin()) {
                        if (l1.above_closed(*b_it)) {
                            w_curr += *w_it;
                        } 
                    } else {
                        assert(angle_it - angles.begin() - 1 >= 0);
                        size_t ix = angle_it - angles.begin() - 1;
                        if (l1.above_closed(*b_it)) {
                            deltas[ix] -= *w_it;
                            w_curr += *w_it;
                        } else {
                            deltas[ix] += *w_it;
                        }
                    }
                    w_it++;
                }
                return w_curr;
            };
            double red_curr = set_delta(red_delta, red.begin(), red.end(), red_w.begin());
            double blue_curr = set_delta(blue_delta, blue.begin(), blue.end(), blue_w.begin());
            for (size_t j = 0; j < angles.size(); j++) {
                double stat = f(red_curr / red_total, blue_curr / blue_total);
                if (max_discrepancy <= stat) {
                    max_line = correct_orientation(p_0, *(j + pb));
                    max_discrepancy = stat;
                }
                red_curr += red_delta[j];
                blue_curr += blue_delta[j];
            }
        }
        return std::make_tuple(max_line, max_discrepancy);
    }


    double unique_accumulate(weight_list const& weights, label_list const& labels) {
        std::vector<bool> label_counts(labels.size(), false);
        double curr = 0;
        for (size_t i = 0; i < labels.size(); i++) {
            if (!label_counts[labels[i]]) {
                curr += weights[i];
                label_counts[labels[i]] = true; 
            }
        }
        return curr;
    }

    std::tuple<Point<>, double> max_halfplane_labeled(
            point_list& point_net,
            point_list& red,
            weight_list& red_w,
            label_list& red_labels,
            point_list& blue,
            weight_list& blue_w,
            label_list& blue_labels,
            std::function<double(double, double)> const& f) {
        /*
        * The maximum label must be less than the maximum number of points. 
        * Labels are assumed to be contigous.
        */
        assert(red.size() == red_w.size());
        assert(blue.size() == blue_w.size());
        assert(blue_labels.size() == blue_w.size());
        assert(red_labels.size() == red_w.size());

        double red_total = unique_accumulate(red_w, red_labels);
        double blue_total = unique_accumulate(blue_w, blue_labels);;

        double max_discrepancy = -std::numeric_limits<double>::infinity();
        Pt2 max_line;

        for (auto p_it = point_net.begin(); p_it != point_net.end() - 1; p_it++) {
            auto p_0 = *p_it;

            auto order_f = [&] (Point<> const& p1) {
                return order_function(p_0, p1);
            };

            auto cmpf = [&](Point<> const& p1, Point<> const& p2) {
                return order_f(p1) < order_f(p2);
            };
            auto pb = p_it + 1;
            std::sort(pb, point_net.end(), cmpf);

            auto l1 = correct_orientation(p_0, *pb);
            std::vector<double> red_delta(point_net.end() - pb, 0.0);
            std::vector<double> blue_delta(point_net.end() - pb, 0.0);

            auto set_delta = [&] (auto& deltas, auto& pts, auto& ws, auto& lbl) {
            
                //Set initial label count and initial weight
                std::vector<int64_t> label_counts(lbl.size(), 0);
                double w_curr = 0;
                {
                    auto bl = lbl.begin(); 
                    auto bp = pts.begin();
                    auto bw = ws.begin();
                    while (bl != lbl.end()) {
                        if (l1.above_closed(*bp)) {
                            if (label_counts[*bl] == 0) {
                                w_curr += *bw;
                            }
                            label_counts[*bl] += 1;
                        } 
                        bl++;
                        bp++;
                        bw++;
                    }
                } 

                // Order these so that they are in order

                std::vector<double> local_orders(pts.size(), 0.0);
                std::transform(pts.begin(), pts.end(), local_orders.begin(), order_f);
                auto permutation = sort_permutation(local_orders, std::less<>());
                //apply_permutation_in_place(orders, permutation);
                apply_permutation_in_place(ws, permutation);
                apply_permutation_in_place(lbl, permutation);
                apply_permutation_in_place(pts, permutation);
                {

                    auto bl = lbl.begin(); 
                    auto bw = ws.begin();
                    auto bp = pts.begin();
                    auto bn = pb;
                    auto bd = deltas.begin();
                    //Set to first line
                    while (order_f(*bp) < order_f(*bn)) {
                        bl++; bw++; bp++;
                    }
                    bn++;
                    //Now walk through the lines.
                    for (; bn != point_net.end(); ++bn, ++bd) {
                        double delta = 0;
                        for (; bp != pts.end() && order_f(*bp) < order_f(*bn); bl++, bw++, bp++) {
                            if (l1.above_closed(*bp)) {
                                //Removing points
                                label_counts[*bl] -= 1;
                                if (label_counts[*bl] <= 0) {
                                    delta -= *bw;
                                }
                            } else {
                                //Adding points
                                if (label_counts[*bl] == 0) {
                                    delta += *bw;
                                }
                                label_counts[*bl] += 1;
                            }
                        }
                        *bd = delta;
                    }
                }
                return w_curr;
            };
            double red_curr = set_delta(red_delta, red, red_w, red_labels);
            double blue_curr = set_delta(blue_delta, blue, blue_w, blue_labels);
            auto pn = pb;
            auto rd = red_delta.begin();
            auto bd = blue_delta.begin();
            for (; pn != point_net.end(); ++pn, ++rd, ++bd) {
                double stat = f(red_curr / red_total, blue_curr / blue_total);
                if (max_discrepancy <= stat) {
                    max_line = correct_orientation(p_0, *pn);
                    max_discrepancy = stat;
                }
                red_curr += *rd;
                blue_curr += *bd;
            }
        }
        return std::make_tuple(max_line, max_discrepancy);
    }


    double under_line(Pt2& line, point_list& points, weight_list& weights) {
        double curr = 0;
        auto pt_it = points.begin();
        auto w_it = weights.begin();
        for (; pt_it != points.end(); ++pt_it) {
            if (line.above_closed(*pt_it)) {
                curr += *w_it;
            }
            ++w_it;
        }
        return curr;
    }

    double evaluate_line(Pt2& line, 
                        point_list& red, 
                        weight_list& red_w,
                        point_list& blue,
                        weight_list& blue_w,
                        std::function<double(double, double)> const& f) {
        double red_total = std::accumulate(red_w.begin(), red_w.end(), 0.0, std::plus<>());
        double blue_total = std::accumulate(blue_w.begin(), blue_w.end(), 0.0, std::plus<>());
        double rf = under_line(line, red, red_w) / red_total;
        double bf = under_line(line, blue, blue_w) / blue_total;
        return f(rf, bf);
    }

    std::tuple<Pt2, double> max_halfplane_simple(
            point_list& point_net,
            point_list& red,
            weight_list& red_w,
            point_list& blue,
            weight_list& blue_w,
            std::function<double(double, double)> const& f) {
        assert(red.size() == red_w.size());
        assert(blue.size() == blue_w.size());
       
        double max_f = -std::numeric_limits<double>::infinity();
        Pt2 max_line;

        for (auto pt_it = point_net.begin(); pt_it != point_net.end() - 1; ++pt_it){
            for (auto pt_it2 = pt_it + 1; pt_it2 != point_net.end(); ++pt_it2) {
                auto line = correct_orientation(*pt_it, *pt_it2);
                double new_max = evaluate_line(line, red, red_w, blue, blue_w, f);
                if (max_f <= new_max) {
                    max_f = new_max;
                    max_line = line;
                }
            }
        }
        return std::make_tuple(max_line, max_f);
    }


    double under_line_labeled(Pt2& line, point_list& points, weight_list& weights, label_list& labels) {
        std::vector<bool> label_counts(labels.size(), false);
        double curr = 0;
        auto pt_it = points.begin();
        auto w_it = weights.begin();
        auto l_it = labels.begin();
        for (; pt_it != points.end(); ++pt_it, ++w_it, ++l_it) {
            if (line.above_closed(*pt_it) && !label_counts[*l_it]) {
                curr += *w_it;
                label_counts[*l_it] = true;
            }
        }
        return curr;
    }

    double evaluate_line_labeled(Pt2& line, 
                        point_list& red, 
                        weight_list& red_w,
                        label_list& red_labels,
                        point_list& blue,
                        weight_list& blue_w,
                        label_list& blue_labels,
                        std::function<double(double, double)> const& f) {

        auto rp = red.begin();
        auto bp = blue.begin();
        double red_total = unique_accumulate(red_w, red_labels);
        double blue_total = unique_accumulate(blue_w, blue_labels);
        double rf = under_line_labeled(line, red, red_w, red_labels) / red_total;
        double bf = under_line_labeled(line, blue, blue_w, blue_labels) / blue_total;
        return f(rf, bf);
    }


    std::tuple<Pt2, double> max_halfplane_simple_labeled(
            point_list& point_net,
            point_list& red,
            weight_list& red_w,
            label_list& red_labels,
            point_list& blue,
            weight_list& blue_w,
            label_list& blue_labels,
            std::function<double(double, double)> const& f) {
       
        double max_f = -std::numeric_limits<double>::infinity();
        Pt2 max_line;

        for (auto pt_it = point_net.begin(); pt_it != point_net.end() - 1; ++pt_it){
            for (auto pt_it2 = pt_it + 1; pt_it2 != point_net.end(); ++pt_it2) {
                auto line = correct_orientation(*pt_it, *pt_it2);
                double new_max = evaluate_line_labeled(line, red, red_w, red_labels, blue, blue_w, blue_labels, f);
                if (max_f <= new_max) {
                    max_f = new_max;
                    max_line = line;
                }
            }
        }
        return std::make_tuple(max_line, max_f);
    }

    Point<3> lift_halfspace(Point<2> const& h, Point<3> const& p) {
        /*
         * Maps a halfspace in 2d to one in 3d using the correct lifting based on the single fixed point in affine.
         */
        return {h[0] * p[2], h[1] * p[2], -p[0] * h[0] - p[1] * h[1] - p[3] * h[2], h[2] * p[2]};
    }

    Point<2> dropPoint(Point<3> const& fixed_point, Point<3> const& p1) {
        /*
         * Does an affine transformation from 3 to 2.
         */
        return {p1[0] * fixed_point[2] - fixed_point[0] * p1[2],
                p1[1] * fixed_point[2] - fixed_point[1] * p1[2],
                p1[3] * fixed_point[2] - fixed_point[3] * p1[2]};
    }

    std::tuple<Pt3, double> max_halfspace(
            point3_list& point_net,
            point3_list& red,
            weight_list& red_w,
            point3_list& blue,
            weight_list& blue_w,
            std::function<double(double, double)> const& f) {
        Pt3 max_halfspace;
        double max_stat = -std::numeric_limits<double>::infinity();
        for (auto pt_it = point_net.begin(); pt_it != point_net.end() - 2; ++pt_it) {
            auto curr_pt = *pt_it;
            auto drop = [&] (Pt3 const& pt) {
                return dropPoint(curr_pt, pt);
            };
            std::vector<Pt2> drop_net(point_net.end() - pt_it - 1, Pt2());
            std::vector<Pt2> drop_red(red.size(), Pt2());
            std::vector<Pt2> drop_blue(blue.size(), Pt2());
            std::transform(pt_it + 1, point_net.end(), drop_net.begin(), drop);
            std::transform(red.begin(), red.end(), drop_red.begin(), drop);
            std::transform(blue.begin(), blue.end(), drop_blue.begin(), drop);
            auto pair_mx = max_halfplane(drop_net, drop_red, red_w, drop_blue, blue_w, f);
            if (std::get<1>(pair_mx) >= max_stat) {
                max_stat = std::get<1>(pair_mx);
                max_halfspace = lift_halfspace(std::get<0>(pair_mx), curr_pt);
            }
        }
        return std::make_tuple(max_halfspace, max_stat);
    }


    std::tuple<Pt3, double> max_halfspace_labeled(
            point3_list& point_net,
            point3_list& red,
            weight_list& red_w,
            label_list& red_labels,
            point3_list& blue,
            weight_list& blue_w,
            label_list& blue_labels,
            std::function<double(double, double)> const& f) {
        Pt3 max_halfspace;
        double max_stat = -std::numeric_limits<double>::infinity();
        for (auto pt_it = point_net.begin(); pt_it != point_net.end() - 2; ++pt_it) {
            auto curr_pt = *pt_it;
            auto drop = [&] (Pt3 const& pt) {
                return dropPoint(curr_pt, pt);
            };
            std::vector<Pt2> drop_net(point_net.end() - pt_it - 1, Pt2());
            std::vector<Pt2> drop_red(red.size(), Pt2());
            std::vector<Pt2> drop_blue(blue.size(), Pt2());
            std::transform(pt_it + 1, point_net.end(), drop_net.begin(), drop);
            std::transform(red.begin(), red.end(), drop_red.begin(), drop);
            std::transform(blue.begin(), blue.end(), drop_blue.begin(), drop);
            auto pair_mx = max_halfplane_labeled(drop_net, drop_red, red_w, red_labels, drop_blue, blue_w, blue_labels, f);
            if (std::get<1>(pair_mx) >= max_stat) {
                max_stat = std::get<1>(pair_mx);
                max_halfspace = lift_halfspace(std::get<0>(pair_mx), curr_pt);
            }
        }
        return std::make_tuple(max_halfspace, max_stat);
    }
}
