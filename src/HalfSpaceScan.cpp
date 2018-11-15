#include <algorithm>
#include <queue>
#include <random>

#include "FunctionApprox.hpp"
#include "Range.hpp"
#include "Segment.hpp"
#include "HalfSpaceScan.hpp"

namespace pyscan{



    inline static pt2_t drop_point(const pt3_t& fixed_point, const pt3_t& p) {
        return {p(0) - fixed_point(0) * p(2) / fixed_point(2),
                p(1) - fixed_point(1) * p(2) / fixed_point(2),
                1.0};
    }

    inline static wpt2_t drop_wpoint(const pt3_t& fixed_point, const wpt3_t& p) {
        return {p.get_weight(),
                p(0) - fixed_point(0) * p(2) / fixed_point(2),
                p(1) - fixed_point(1) * p(2) / fixed_point(2),
                1.0};
    }

    inline static lpt2_t drop_lpoint(const pt3_t& fixed_point, const lpt3_t& p) {
        return {p.get_label(),
                p.get_weight(),
                p(0) - fixed_point(0) * p(2) / fixed_point(2),
                p(1) - fixed_point(1) * p(2) / fixed_point(2),
                1.0};
    }

    inline static pt3_t lift_half_space(const pt2_t& h, const pt3_t& p) {
        return {h[0] * p[2],
                h[1] * p[2],
                -p[0] * h[0] - p[1] * h[1] - p[3] * h[2],
                h[2] * p[2]};
    }

    inline double plane_order(Point<2> const& p1, Point<2> const& p2) {
        double a = util::det2(p1[1], p1[2], p2[1], p2[2]);
        double b = -util::det2(p1[0], p1[2], p2[0], p2[2]);
        double orientation = -std::copysign(1.0, b);
        double inv_norm = 1 / sqrt(a * a + b * b);
        return -a * inv_norm * orientation;
    }

    std::tuple<halfspace2_t, double> max_halfplane(
            const point_list_t& point_net,
            const wpoint_list_t& red,
            const wpoint_list_t& blue,
            const discrepancy_func_t& f) {

        double red_total = computeTotal(red);
        double blue_total = computeTotal(blue);

        double max_discrepancy = -std::numeric_limits<double>::infinity();
        halfspace2_t max_plane;

        assert(point_net.size() >= 2);
        for (size_t i = 0; i < point_net.size() - 1; ++i) {
            auto pivot = point_net[i];

            std::vector<halfspace2_t> halfplanes;
            for (size_t j = i + 1; j < point_net.size(); ++j) {
                halfplanes.emplace_back(pivot, point_net[j]);
            }
            std::sort(halfplanes.begin(), halfplanes.end(), [](const halfspace2_t& p1, const halfspace2_t& p2) {
                return -p1[0] < -p2[0];
            });

            auto& l1 = halfplanes[0];
            std::vector<double> red_delta(halfplanes.size() - 1, 0.0);
            std::vector<double> blue_delta(halfplanes.size() - 1, 0.0);
            std::vector<double> angles;
            angles.reserve(halfplanes.size());
            for (auto& plane : halfplanes) {
                angles.emplace_back(-plane[0]);
            }
            auto calc_delta = [&](const wpoint_list_t& pts, std::vector<double>& deltas) {
                double res = 0.0;
                for (auto const& pt : pts) {
                    double val = plane_order(pivot, pt);
                    auto angle_it = std::lower_bound(angles.begin(), angles.end(), val);
                    //If the angle is begin or end then it is in the last wedge and we don't count it.
                    if (angle_it == angles.end() || angle_it == angles.begin()){
                        if (l1.contains(pt)) {
                            res += pt.get_weight();
                        }
                    } else {
                        auto ix = std::distance(angles.begin(), angle_it) - 1;
                        if (l1.contains(pt)) {
                            deltas[ix] -= pt.get_weight();
                            res += pt.get_weight();
                        } else {
                            deltas[ix] += pt.get_weight();
                        }
                    }
                }
                return res;
            };

            double red_curr = calc_delta(red, red_delta);
            double blue_curr = calc_delta(blue, blue_delta);
            for (size_t j = 0; true ;++j) {
                double stat = f(red_curr, red_total, blue_curr, blue_total);
                if (max_discrepancy <= stat) {
                    max_plane = halfplanes[j];
                    max_discrepancy = stat;
                }
                if (j == halfplanes.size() - 1) {
                    break;
                }
                red_curr += red_delta[j];
                blue_curr += blue_delta[j];
            }
        }

        return std::make_tuple(max_plane, max_discrepancy);
    }


    std::tuple<halfspace3_t, double> max_halfspace(
            const point3_list_t& point_net,
            const wpoint3_list_t& red,
            const wpoint3_list_t& blue,
            const discrepancy_func_t& f) {

        pt3_t max_halfspace;
        double max_discrepancy = -std::numeric_limits<double>::infinity();

        assert(point_net.size() >= 3);
        for (size_t i = 0; i < point_net.size() - 2; ++i) {
            auto pivot = point_net[i];
            auto drop = [&pivot](const pt3_t& pt) {
                return drop_point(pivot, pt);
            };
            auto wdrop = [&pivot](const wpt3_t& pt) {
                return drop_wpoint(pivot, pt);
            };
            point_list_t drop_net(point_net.size() - i - 1, pt2_t());
            wpoint_list_t drop_red(red.size(), wpt2_t());
            wpoint_list_t drop_blue(blue.size(), wpt2_t());
            std::transform(point_net.begin() + i + 1, point_net.end(), drop_net.begin(), drop);
            std::transform(red.begin(), red.end(), drop_red.begin(), wdrop);
            std::transform(blue.begin(), blue.end(), drop_blue.begin(), wdrop);
            auto pair_mx = max_halfplane(drop_net, drop_red, drop_blue, f);
            if (std::get<1>(pair_mx) >= max_discrepancy) {
                max_discrepancy = std::get<1>(pair_mx);
                max_halfspace = lift_half_space(std::get<0>(pair_mx).get_coords(), pivot);
            }
        }
        return std::make_tuple(halfspace3_t(max_halfspace), max_discrepancy);
    }

    struct LabeledValue {
        size_t label;
        double value;
    };

    using wedge_t = std::vector<LabeledValue>;
    using label_set_t = std::unordered_map<size_t , size_t >;

    inline static double update_weight(
            std::unordered_map<size_t, size_t> &cur_set,
            const wedge_t &adding, const wedge_t& removing) {

        double update_diff = 0.0;
        for (auto &x: adding) {
            auto it = cur_set.find(x.label);
            if (it != cur_set.end()) {
                cur_set[x.label]++;
            } else {
                update_diff += x.value;
                cur_set[x.label] = 1;
            }
        }

        for (auto &x: removing) {
            auto it = cur_set.find(x.label);
            assert(it != cur_set.end());
            assert(it->second > 0);

            if (it->second == 1) {
                update_diff -= x.value;
                cur_set.erase(it);
            } else {
                cur_set[x.label]--;
            }
        }

        return update_diff;
    }

    std::tuple<halfspace2_t, double> max_halfplane_labeled(
            const point_list_t& point_net,
            const lpoint_list_t& red,
            const lpoint_list_t& blue,
            const discrepancy_func_t& f) {

        double red_total = computeTotal(red);
        double blue_total = computeTotal(blue);

        double max_discrepancy = -std::numeric_limits<double>::infinity();
        halfspace2_t max_plane;

        assert(point_net.size() >= 2);
        for (size_t i = 0; i < point_net.size() - 1; ++i) {
            auto pivot = point_net[i];

            std::vector<halfspace2_t> halfplanes;
            for (size_t j = i + 1; j < point_net.size(); ++j) {
                halfplanes.emplace_back(pivot, point_net[j]);
            }
            std::sort(halfplanes.begin(), halfplanes.end(), [](const halfspace2_t& p1, const halfspace2_t& p2) {
                return -p1[0] < -p2[0];
            });

            auto& l1 = halfplanes[0];

            std::vector<double> angles;
            angles.reserve(halfplanes.size());
            for (auto& plane : halfplanes) {
                angles.emplace_back(-plane[0]);
            }

            std::vector<wedge_t> red_deltaR(halfplanes.size() - 1), blue_deltaR(halfplanes.size() - 1);
            std::vector<wedge_t> red_deltaA(halfplanes.size() - 1), blue_deltaA(halfplanes.size() - 1);

            auto calc_delta = [&](const lpoint_list_t& pts,
                    std::vector<wedge_t>& deltaR,
                    std::vector<wedge_t>& deltaA,
                    label_set_t& labels) {
                double res = 0.0;
                for (auto const& pt : pts) {

                    auto angle_it = std::lower_bound(angles.begin(), angles.end(), plane_order(pivot, pt));
                    //If the angle is begin or end then it is in the last wedge and we don't count it.
                    if (l1.contains(pt)) {
                        auto it = labels.find(pt.get_label());
                        if (it == labels.end()) {
                            res += pt.get_weight();
                            labels[pt.get_label()] = 1;
                        } else {
                            labels[pt.get_label()]++;
                        }
                    }
                    if (angle_it == angles.end() || angle_it == angles.begin()){
                        continue;
                    } else {
                        auto ix = std::distance(angles.begin(), angle_it) - 1;
                        if (l1.contains(pt)) {
                            deltaR[ix].emplace_back(LabeledValue{pt.get_label(), pt.get_weight()});
                        } else {
                            deltaA[ix].emplace_back(LabeledValue{pt.get_label(), pt.get_weight()});
                        }
                    }
                }
                return res;
            };

            label_set_t red_set;
            label_set_t blue_set;
            double red_curr = calc_delta(red, red_deltaR, red_deltaA, red_set);
            double blue_curr = calc_delta(blue, blue_deltaR, blue_deltaA, blue_set);
            for (size_t j = 0; true ;++j) {
                double new_stat = f(red_curr, red_total, blue_curr, blue_total);
                if (max_discrepancy <= new_stat) {
                    max_plane = halfplanes[j];
                    max_discrepancy = new_stat;
                }
                if (j == halfplanes.size() - 1) {
                    break;
                }
                red_curr += update_weight(red_set, red_deltaA[j], red_deltaR[j]);
                blue_curr += update_weight(blue_set, blue_deltaA[j], blue_deltaR[j]);
            }
        }

        return std::make_tuple(max_plane, max_discrepancy);
    }


    template <typename T>
    void remove_duplicates_alt(T& pts) {
        std::sort(pts.begin(), pts.end(), [](auto const& p1, auto const& p2){
            return p1[0] < p2[1];

        });

        auto end_it = std::unique(pts.begin(), pts.end(), [] (auto const& p1, auto const& p2) {
            return util::aeq(std::abs(p1[0] - p2[0]) + std::abs(p1[1] - p2[1]), 0.0);
        });
        pts.erase(end_it, pts.end());
    }

     void approx_hull(lpoint_it_t begin, lpoint_it_t end, double eps, lpoint_list_t& output) {
        if (begin == end) {
            return ;
        }
        auto first_pt = *begin;
        auto max_f = [&] (Vec2 direction) {
            double max_dir = -std::numeric_limits<double>::infinity();
            Vec2 curr_pt {0.0, 0.0};
            for (auto b = begin;  b != end; b++) {
                auto& pt = *b;
                double curr_dir = direction[0] * pt(0) + direction[1] * pt(1);
                if (max_dir < curr_dir) {
                    max_dir = curr_dir;
                    curr_pt[0] = pt(0);
                    curr_pt[1] = pt(1);
                }
            }
            return curr_pt;
        };
        auto vecs = eps_core_set(eps, max_f);
        remove_duplicates_alt(vecs);
        if (vecs.size() > static_cast<size_t>(end - begin)) {
            for (; begin != end; begin++) {
                output.emplace_back(*begin);
            }
        } else {
            for (auto &v: vecs) {
                //Assume that all of the begin to end has same label and weight.
                output.emplace_back(LPoint(first_pt.get_label(), first_pt.get_weight(), v[0], v[1], 1.0));
            }
        }
    }

    lpoint_list_t filter_kernel2d(lpoint_list_t lpoints, double alpha) {


        auto comp = [] (lpt2_t const& p1, lpt2_t const& p2) {
            return p1.get_label() < p2.get_label();
        };
        std::sort(lpoints.begin(), lpoints.end(), comp) ;

        lpoint_list_t output_pts(lpoints.size(), lpt2_t());
        for (auto pt_it = lpoints.begin(); pt_it != lpoints.end(); ) {
            auto end_of_group = std::upper_bound(pt_it, lpoints.end(), *pt_it, comp);
            approx_hull(pt_it, end_of_group, alpha, output_pts);
        }
        return output_pts;

    }

    std::tuple<halfspace2_t, double> max_halfplane_labeled_restricted(
            const point_list_t& point_net,
            const lpoint_list_t& red,
            const lpoint_list_t& blue,
            double red_total,
            double blue_total,
            const filter_func2_t& filter,
            const discrepancy_func_t& f) {

        double max_discrepancy = 0;
        halfspace2_t max_plane;

        assert(point_net.size() >= 2);
        for (size_t i = 0; i < point_net.size() - 1; ++i) {
            auto pivot = point_net[i];

            std::vector<halfspace2_t> halfplanes;
            for (size_t j = i + 1; j < point_net.size(); ++j) {
                halfspace2_t possible(pivot, point_net[j]);
                if (filter(possible)) {
                    halfplanes.emplace_back(possible);
                }
            }
            if (halfplanes.empty()) {
                continue;
            }

            std::sort(halfplanes.begin(), halfplanes.end(), [](const halfspace2_t& p1, const halfspace2_t& p2) {
                return -p1[0] < -p2[0];
            });

            auto& l1 = halfplanes[0];

            std::vector<double> angles;
            angles.reserve(halfplanes.size());
            for (auto& plane : halfplanes) {
                angles.emplace_back(-plane[0]);
            }

            std::vector<wedge_t> red_deltaR(halfplanes.size() - 1), blue_deltaR(halfplanes.size() - 1);
            std::vector<wedge_t> red_deltaA(halfplanes.size() - 1), blue_deltaA(halfplanes.size() - 1);

            auto calc_delta = [&](const lpoint_list_t& pts,
                                  std::vector<wedge_t>& deltaR,
                                  std::vector<wedge_t>& deltaA,
                                  label_set_t& labels) {
                double res = 0.0;
                for (auto const& pt : pts) {

                    auto angle_it = std::lower_bound(angles.begin(), angles.end(), plane_order(pivot, pt));
                    //If the angle is begin or end then it is in the last wedge and we don't count it.
                    if (l1.contains(pt)) {
                        auto it = labels.find(pt.get_label());
                        if (it == labels.end()) {
                            res += pt.get_weight();
                            labels[pt.get_label()] = 1;
                        } else {
                            labels[pt.get_label()]++;
                        }
                    }
                    if (angle_it == angles.end() || angle_it == angles.begin()){
                        continue;
                    } else {
                        auto ix = std::distance(angles.begin(), angle_it) - 1;
                        if (l1.contains(pt)) {
                            deltaR[ix].emplace_back(LabeledValue{pt.get_label(), pt.get_weight()});
                        } else {
                            deltaA[ix].emplace_back(LabeledValue{pt.get_label(), pt.get_weight()});
                        }
                    }
                }
                return res;
            };

            label_set_t red_set;
            label_set_t blue_set;
            double red_curr = calc_delta(red, red_deltaR, red_deltaA, red_set);
            double blue_curr = calc_delta(blue, blue_deltaR, blue_deltaA, blue_set);
            for (size_t j = 0; true ;++j) {
                double new_stat = f(red_curr, red_total, blue_curr, blue_total);
                if (max_discrepancy <= new_stat) {
                    max_plane = halfplanes[j];
                    max_discrepancy = new_stat;
                }
                if (j == halfplanes.size() - 1) {
                    break;
                }
                red_curr += update_weight(red_set, red_deltaA[j], red_deltaR[j]);
                blue_curr += update_weight(blue_set, blue_deltaA[j], blue_deltaR[j]);
            }
        }

        return std::make_tuple(max_plane, max_discrepancy);
    }

    std::tuple<halfspace3_t, double> max_halfspace_labeled(
            const point3_list_t& point_net,
            const lpoint3_list_t& red,
            const lpoint3_list_t& blue,
            const discrepancy_func_t& f) {
        pt3_t max_halfspace;
        double max_discrepancy = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < point_net.size() - 2; ++i) {
            auto pivot = point_net[i];
            auto drop = [&pivot](const pt3_t& pt) {
                return drop_point(pivot, pt);
            };
            auto ldrop = [&pivot](const lpt3_t& pt) {
                return drop_lpoint(pivot, pt);
            };
            point_list_t drop_net(point_net.size() - i - 1, pt2_t());
            lpoint_list_t drop_red(red.size(), lpt2_t());
            lpoint_list_t drop_blue(blue.size(), lpt2_t());
            std::transform(point_net.begin() + i + 1, point_net.end(), drop_net.begin(), drop);
            std::transform(red.begin(), red.end(), drop_red.begin(), ldrop);
            std::transform(blue.begin(), blue.end(), drop_blue.begin(), ldrop);
            auto pair_mx = max_halfplane_labeled(drop_net, drop_red, drop_blue, f);
            if (std::get<1>(pair_mx) >= max_discrepancy) {
                max_discrepancy = std::get<1>(pair_mx);
                max_halfspace = lift_half_space(std::get<0>(pair_mx).get_coords(), pivot);
            }
        }
        return std::make_tuple(halfspace3_t(max_halfspace), max_discrepancy);
    }

    std::tuple<halfspace3_t, double> max_halfspace_labeled_restricted(
            const pt3_t& pt,
            const lpoint3_list_t& point_net,
            const lpoint3_list_t& red,
            const lpoint3_list_t& blue,
            double red_total, double blue_total,
            double alpha,
            const filter_func3_t& filter,
            const discrepancy_func_t& f) {

        auto pivot = pt;
        auto drop = [&pivot](const pt3_t& pt) {
            return drop_point(pivot, pt);
        };
        auto ldrop = [&pivot](const lpt3_t& pt) {
            return drop_lpoint(pivot, pt);
        };
        point_list_t drop_net(point_net.size(), pt2_t());
        lpoint_list_t drop_red(red.size(), lpt2_t());
        lpoint_list_t drop_blue(blue.size(), lpt2_t());
        std::transform(point_net.begin(), point_net.end(), drop_net.begin(), drop);
        std::transform(red.begin(), red.end(), drop_red.begin(), ldrop);
        std::transform(blue.begin(), blue.end(), drop_blue.begin(), ldrop);

        auto func = [filter, pivot] (halfspace2_t const& h) {
            HalfSpace<3> halfspace_tmp(lift_half_space(h.get_coords(), pivot));
            return filter(halfspace_tmp);
        };

        auto [h, hmx] = max_halfplane_labeled_restricted(drop_net, drop_red, drop_blue, red_total, blue_total, func, f);
        return std::make_tuple(halfspace3_t(lift_half_space(h.get_coords(), pivot)), hmx);
    }


    std::tuple<halfspace2_t, double> max_halfplane_simple(
            const point_list_t &point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            const discrepancy_func_t &f) {
        return max_range2<halfspace2_t, WPoint>(point_net, red, blue, f);
    }


    std::tuple<halfspace3_t, double> max_halfspace_simple(
            const point3_list_t &point_net,
            const wpoint3_list_t &red,
            const wpoint3_list_t &blue,
            const discrepancy_func_t &f) {
        return max_range3<halfspace3_t, WPoint, 3>(point_net, red, blue, f);
    }

    std::tuple<halfspace2_t, double> max_halfplane_labeled_simple(
            const point_list_t &point_net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            const discrepancy_func_t &f) {
        return max_range2<halfspace2_t, LPoint>(point_net, red, blue, f);
    }

    std::tuple<halfspace3_t, double> max_halfspace_labeled_simple(
            const point3_list_t &point_net,
            const lpoint3_list_t &red,
            const lpoint3_list_t &blue,
            const discrepancy_func_t &f) {
        return max_range3<halfspace3_t, LPoint, 3>(point_net, red, blue, f);
    }


    using seg_t = Segment;
    using seg_list_t = std::vector<Segment>;

    class CutCell {
        double red_w = 0;
        double blue_w = 0;
        double cell_sc = 1.0;
        seg_list_t red_segments;
        seg_list_t blue_segments;
    public:

        CutCell() {}

        CutCell(const point_list_t& red, const point_list_t& blue) {
            for (auto const& pt : red) {
                red_segments.emplace_back(pt);
            }
            for (auto const& pt : blue) {
                blue_segments.emplace_back(pt);
            }
        }

        CutCell(seg_list_t red, double rw, seg_list_t blue, double bw, double csc) :
            red_w(rw), blue_w(bw), cell_sc(csc),
            red_segments(std::move(red)), blue_segments(std::move(blue)) {}

        double red_weight() const { return red_w; }
        double blue_weight() const { return red_w; }
        double cell_scale() const { return cell_sc; }

        template <typename RNG>
        std::tuple<CutCell, CutCell> split(RNG& rng) const {
            // Uniformly

            Segment seg;
            std::uniform_real_distribution<double> dist(0, 1);
            if ((blue_segments.empty() || dist(rng) < .5) && !red_segments.empty()) {
                std::uniform_int_distribution<size_t> distribution(0, red_segments.size() - 1);
                seg = red_segments[distribution(rng)];
            } else if (!blue_segments.empty()){
                std::uniform_int_distribution<size_t> distribution(0, blue_segments.size() - 1);
                seg = blue_segments[distribution(rng)];
            } else {
                //Return two empty cells.
                return {CutCell(), CutCell()};
            }
            auto split_segments = [&](seg_list_t const& crossing_segments) {

                double below = 0;
                seg_list_t upper_crossing_seg;
                seg_list_t lower_crossing_seg;

                for (auto& tseg : crossing_segments) {
                    if (tseg.approx_eq(seg)) {
                        below += cell_sc;
                        continue;
                    } else if (tseg.crossed(seg)) {
                        auto [upper, lower] = tseg.split(seg);
                        lower_crossing_seg.emplace_back(lower);
                        upper_crossing_seg.emplace_back(upper);
                    } else if (tseg.gte(seg)) {
                        upper_crossing_seg.push_back(tseg);
                    } else {
                        lower_crossing_seg.push_back(tseg);
                        below += cell_sc;
                    }
                }
                return make_tuple(upper_crossing_seg, lower_crossing_seg, below);
            };

            auto [u_red_crossing, l_red_crossing, red_below_w] = split_segments(red_segments);
            auto [u_blue_crossing, l_blue_crossing, blue_below_w] = split_segments(blue_segments);

            //Subsample the lower and upper lines.
//            std::shuffle(u_red_crossing.begin(), u_red_crossing.end(), rng);
//            std::shuffle(l_red_crossing.begin(), l_red_crossing.end(), rng);
//            std::shuffle(u_blue_crossing.begin(), u_blue_crossing.end(), rng);
//            std::shuffle(l_blue_crossing.begin(), l_blue_crossing.end(), rng);
//            size_t new_cell_size = u_red_crossing.size() + l_red_crossing.size() +
//                    u_blue_crossing.size() + l_blue_crossing.size();
//
//            //double u_prop = (u_red_crossing.size() + u_blue_crossing.size()) / static_cast<double>(new_cell_size);
//            //double l_prop = 1 - u_prop;
//
//            double rescale = static_cast<double>(red_segments.size() + blue_segments.size()) / new_cell_size;
//
//            auto resize = [rescale](seg_list_t& pts) {
//                size_t new_size = lround(rescale * pts.size());
//                if (new_size >= pts.size()) {
//                    return ;
//                } else {
//                    pts.resize(new_size);
//                }
//            };
//
//            size_t ursz = u_red_crossing.size(), ubsz = u_blue_crossing.size(),
//                lrsz = l_red_crossing.size(), lbsz = l_blue_crossing.size();
//
//            resize(u_red_crossing);
//            resize(l_red_crossing);
//            resize(u_blue_crossing);
//            resize(l_blue_crossing);
//
//            double upper_scale = (ursz + ubsz) / static_cast<double>(u_red_crossing.size() + u_blue_crossing.size());
//            double lower_scale = (lrsz + lbsz) / static_cast<double>(l_red_crossing.size() + l_blue_crossing.size());
            //Need the red lines and the blue lines to represent the original sets.
//            CutCell c1(u_red_crossing, red_below_w + red_w,
//                        u_blue_crossing, blue_below_w + blue_w,
//                        cell_sc * upper_scale);
//            CutCell c2(l_red_crossing, red_w, l_blue_crossing, blue_w, cell_sc * lower_scale);
            CutCell c1(u_red_crossing, red_w + red_below_w,
                       u_blue_crossing, blue_w + blue_below_w,
                       1.0);
            CutCell c2(l_red_crossing, red_w, l_blue_crossing, blue_w, 1.0);
            return {c1, c2};
        }

        template<typename RNG>
        std::tuple<Point<2>, double, double> choose_line(RNG& rng) const {
            Segment seg;
            std::uniform_real_distribution<double> dist(0, 1);
            if ((blue_segments.empty() || dist(rng) < .5) && !red_segments.empty()) {
                std::uniform_int_distribution<size_t> distribution(0, red_segments.size() - 1);
                seg = red_segments[distribution(rng)];
            } else if (!blue_segments.empty()){
                std::uniform_int_distribution<size_t> distribution(0, blue_segments.size() - 1);
                seg = blue_segments[distribution(rng)];
            } else {
                //Return the line at infinity
                return {Point<2>(0.0, 0.0, 1.0), 0.0, 0.0};
            }
            //Choose an endpoint
            auto line = seg.get_e1().orient_down();

            double red_l_weight = red_w;
            for (auto& pt : red_segments) {
                if (line.above_closed(pt)) red_l_weight += cell_sc;
            }

            double blue_l_weight = blue_w;
            for (auto& pt : blue_segments) {
                if (line.above_closed(pt)) blue_l_weight += cell_sc;
            }
            return {line, red_l_weight, blue_l_weight};
        }

        double get_total_weight() const {
            return cell_sc * (red_segments.size() + blue_segments.size());
        }

        bool operator<(const CutCell& s2) const {
            return get_total_weight() < s2.get_total_weight();
        }

        seg_list_t get_red_set() const {
            return red_segments;
        }

        seg_list_t get_blue_set() const {
            return blue_segments;
        }
    };

    std::tuple<halfspace2_t, double> max_halfplane_fast(size_t plane_count,
                                                        const point_list_t &red,
                                                        const point_list_t &blue,
                                                        const discrepancy_func_t &f) {


        auto red_total = static_cast<double>(red.size());
        auto blue_total = static_cast<double>(blue.size());
        std::random_device rd;
        std::minstd_rand gen(rd());
        std::priority_queue<CutCell> cells;

        cells.emplace(red, blue);

        while (cells.size() < plane_count) {

            auto& curr_node = cells.top();

            if (curr_node.get_total_weight() < 2) {
                break;
            }
            auto [upper_cell, lower_cell] = curr_node.split(gen);
            cells.pop();
            cells.push(upper_cell);
            cells.push(lower_cell);
        }


        //Now just choose a single vertex from each cell and compute its discrepancy.
        double max_disc = 0;
        HalfSpace<2> max_halfspace;
        while (!cells.empty()) {
            auto& curr_node = cells.top();
            auto [line, red_w, blue_w] = curr_node.choose_line(gen);
            double scan_val = f(red_w, red_total, blue_w, blue_total);
            if (scan_val > max_disc) {
                max_halfspace = HalfSpace<2>(line);
                max_disc = scan_val;
            }
            cells.pop();
        }
        return {max_halfspace, max_disc};
    }
}
