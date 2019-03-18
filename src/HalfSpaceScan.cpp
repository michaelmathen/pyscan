#include <algorithm>
#include <queue>
#include <random>

#include "FunctionApprox.hpp"
#include "Range.hpp"
#include "Segment.hpp"
#include "ConvexHull.hpp"
#include "HalfSpaceScan.hpp"

namespace pyscan{

    inline static pt2_t internal_drop_point(const pt3_t& fixed_point, const pt3_t& p) {
        return pt2_t(fixed_point[0] * p[1] - fixed_point[1] * p[0],
                     fixed_point[0] * p[2] - fixed_point[2] * p[0],
                     fixed_point[0] * p[3] - fixed_point[3] * p[0]
        );
    }


    pt2_t drop_point(const pt3_t& fixed_point, const pt3_t& p) {
        return internal_drop_point(fixed_point, p);
    }

    wpt2_t drop_point(const pt3_t& fixed_point, const wpt3_t& p) {
       auto new_coords = internal_drop_point(fixed_point, p);
       return wpt2_t(p.get_weight(), new_coords[0], new_coords[1], new_coords[2]);
    }

    lpt2_t drop_point(const pt3_t& fixed_point, const lpt3_t& p) {
        auto new_coords = internal_drop_point(fixed_point, p);
        return lpt2_t(p.get_label(), p.get_weight(), new_coords[0], new_coords[1], new_coords[2]);
    }

//    template <typename Pt>
//    point_list_t drop_points_lpoints(const pt3_t& fixed_pt, const std::vector<Pt>& pts, double& x_offset, ) {
//        x_offset = 1.0 - fixed_pt[0];
//        for (auto& pt : pts) {
//            if ()
//            Pt p_tmp = pt;
//            p_tmp[0] += x_offset;
//
//        }
//
//
//    }

    pt3_t lift_half_space(const halfspace2_t& h, const pt3_t& p) {
        return pt3_t(-(h[0] * p[1] + h[1] * p[2] + h[2] * p[3]) ,
                h[0] * p[0],
                h[1] * p[0],
                h[2] * p[0] );
    }

    inline double plane_order(Point<2> const& p1, Point<2> const& p2) {
        double a = util::det2(p1[1], p1[2], p2[1], p2[2]);
        double b = -util::det2(p1[0], p1[2], p2[0], p2[2]);
        double orientation = -std::copysign(1.0, b);
        double inv_norm = 1 / sqrt(a * a + b * b);
        return -a * inv_norm * orientation;
    }

    std::tuple<halfspace2_t, double> max_halfplane_internal(
            const point_list_t& point_net,
            const wpoint_list_t& red,
            const wpoint_list_t& blue,
            const filter_func2_t& filter,
            const discrepancy_func2_t& f) {

        double max_discrepancy = -std::numeric_limits<double>::infinity();
        halfspace2_t max_plane;

        assert(point_net.size() >= 2);
        for (size_t i = 0; i < point_net.size() - 1; ++i) {
            auto pivot = point_net[i];

            std::vector<halfspace2_t> halfplanes;
            for (size_t j = i + 1; j < point_net.size(); ++j) {
                if (!pivot.approx_eq(point_net[j])) {
                    halfspace2_t possible(pivot, point_net[j]);
                    if (filter(possible) && !std::isnan(possible[0]) && !std::isnan(possible[1]) &&
                        !std::isnan(possible[2])) {
                        halfplanes.push_back(possible);
                    }
                }
            }
            std::sort(halfplanes.begin(), halfplanes.end(), [](const halfspace2_t& p1, const halfspace2_t& p2) {
                return -p1[0] < -p2[0];
            });

            auto new_end = std::unique(halfplanes.begin(), halfplanes.end(), [](const halfspace2_t& p1, const halfspace2_t& p2){
                return util::aeq(p1[0], p2[0]);
            });
            halfplanes.resize( std::distance(halfplanes.begin(), new_end));

            if (halfplanes.empty()) {
                continue;
            }

            auto l1 = halfplanes[0];
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
                    if (l1.contains(pt)) {
                        res += pt.get_weight();
                    }

                    halfspace2_t h_tmp(pivot, pt);
                    if (!std::isnan(h_tmp[0])) {
                        auto angle_it = std::lower_bound(angles.begin(), angles.end(), -h_tmp[0]);
                        //If the angle is begin or end then it is in the last wedge and we don't count it.
                        if (!(angle_it == angles.end() || angle_it == angles.begin())) {
                            auto ix = std::distance(angles.begin(), angle_it) - 1;
                            if (l1.contains(pt)) {
                                deltas[ix] -= pt.get_weight();
                            } else {
                                deltas[ix] += pt.get_weight();
                            }
                        }
                    }
                }
                return res;
            };

            double red_curr = calc_delta(red, red_delta);
            double blue_curr = calc_delta(blue, blue_delta);
            for (size_t j = 0; true ;++j) {
                double stat = f(red_curr, blue_curr);
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


    struct LabeledValue {
        size_t label;
        double value;
    };

    std::ostream &operator<<(std::ostream& os, LabeledValue& val) {
        os << "L(" << val.label << ", "  << val.value << ")";
        return os;
    }

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


    std::tuple<halfspace2_t, double> max_halfplane_internal(
            const point_list_t& point_net,
            const lpoint_list_t& red,
            const lpoint_list_t& blue,
            const filter_func2_t& filter,
            const discrepancy_func2_t& f) {

        double max_discrepancy = 0.0;
        halfspace2_t max_plane;

        if (point_net.size() < 2) {
            return {max_plane, max_discrepancy};
        }
        for (size_t i = 0; i < point_net.size() - 1; ++i) {
            auto pivot = point_net[i];

            std::vector<halfspace2_t> halfplanes;
            for (size_t j = i + 1; j < point_net.size(); ++j) {
                if (!pivot.approx_eq(point_net[j])) {
                    halfspace2_t possible(pivot, point_net[j]);
                    if (filter(possible) && !std::isnan(possible[0]) && !std::isnan(possible[1]) &&
                        !std::isnan(possible[2])) {
                        halfplanes.push_back(possible);
                    }
                }
            }
            std::sort(halfplanes.begin(), halfplanes.end(), [](const halfspace2_t& p1, const halfspace2_t& p2) {
                return -p1[0] < -p2[0];
            });
            auto new_end = std::unique(halfplanes.begin(), halfplanes.end(), [](const halfspace2_t& p1, const halfspace2_t& p2){
                return util::aeq(p1[0], p2[0]);
            });
            halfplanes.resize( std::distance(halfplanes.begin(), new_end));

            if (halfplanes.empty()) {
                continue;
            }

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

                    halfspace2_t h_tmp(pivot, pt);
                    if (!std::isnan(h_tmp[0])) {
                        auto angle_it = std::lower_bound(angles.begin(), angles.end(), -h_tmp[0]);

                        if (angle_it == angles.end() || angle_it == angles.begin()) {
                            continue;
                        } else {
                            auto ix = std::distance(angles.begin(), angle_it) - 1;
                            if (l1.contains(pt)) {
                                deltaR[ix].push_back(LabeledValue{pt.get_label(), pt.get_weight()});
                            } else {
                                deltaA[ix].push_back(LabeledValue{pt.get_label(), pt.get_weight()});
                            }
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
                double new_stat = f(red_curr, blue_curr);
//                if (!util::aeq(range_weight(halfplanes[j], red), red_curr)) {
//                    std::cout << halfplanes[j].str() << std::endl;
//                    std::cout << range_weight(halfplanes[j], red) << " " <<   red_curr << std::endl;
//                    std::cout << update_weight(red_set, red_deltaA[j], red_deltaR[j]) << std::endl;
//                    std::cout << red_deltaA[j - 1] << std::endl;
//                    std::cout << red_deltaR[j - 1] << std::endl;
//                    std::cout << angles << std::endl;
//                    for (auto& pt : red) {
//                        if (halfplanes[j].contains(pt) != halfplanes[j - 1].contains(pt) ) {
//                            std::cout << pt << std::endl;
//                            auto h_tmp = halfspace2_t(pivot, pt);
//                            std::cout << h_tmp.str() << std::endl;
//                            std::cout << *(std::lower_bound(angles.begin(), angles.end(), -h_tmp[0])) << std::endl;
//                        }
//                    }
//                    assert(util::aeq(range_weight(halfplanes[j], red), red_curr));
//                    assert(util::aeq(range_weight(halfplanes[j], blue), blue_curr));
//                }


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


    std::tuple<halfspace2_t, double> max_halfplane(
            const point_list_t& point_net,
            const wpoint_list_t& red,
            const wpoint_list_t& blue,
            const discrepancy_func_t& f) {
        double m_Total = computeTotal(red);
        double b_Total = computeTotal(blue);

        return max_halfplane_internal(point_net, red, blue,
                                      [](halfspace2_t const& h) { (void)h; return true;},
                                      [&](double m, double b) { return f(m, m_Total, b, b_Total);}
        );
    }

    std::tuple<halfspace2_t, double> max_halfplane_labeled(
            const point_list_t& point_net,
            const lpoint_list_t& red,
            const lpoint_list_t& blue,
            const discrepancy_func_t& f) {
        double m_Total = computeTotal(red);
        double b_Total = computeTotal(blue);

        return max_halfplane_internal(point_net, red, blue,
                                      [](halfspace2_t const& h) { (void)h; return true;},
                                      [&](double m, double b) { return f(m, m_Total, b, b_Total);}
        );
    }


    template<template <int> typename P=WPoint>
    std::tuple<halfspace3_t, double> max_halfspace_internal(
            const point3_list_t& point_net,
            const std::vector<P<3>>& red,
            const std::vector<P<3>>& blue,
            const discrepancy_func_t& f) {

        double m_Total = computeTotal(red);
        double b_Total = computeTotal(blue);

        HalfSpace<3> max_halfspace;
        double max_discrepancy = 0.0;

        if (point_net.size() < 3) {
            return {max_halfspace, 0.0};
        }
        for (size_t i = 0; i < point_net.size() - 2; ++i) {
            auto pivot = point_net[i];
            auto drop = [&pivot](const pt3_t& pt) {
                return drop_point(pivot, pt);
            };
            auto wdrop = [&pivot](const P<3>& pt) {
                assert(!pivot.approx_eq(pt));
                return drop_point(pivot, pt);
            };
            point_list_t drop_net(point_net.size() - i - 1, pt2_t());
            std::vector<P<2>> drop_red(red.size(), P<2>());
            std::vector<P<2>> drop_blue(blue.size(), P<2>());
            std::transform(point_net.begin() + i + 1, point_net.end(), drop_net.begin(), drop);
            std::transform(red.begin(), red.end(), drop_red.begin(), wdrop);
            std::transform(blue.begin(), blue.end(), drop_blue.begin(), wdrop);
            //std::cout << blue[0] << drop_blue[0] << std::endl;
            auto [h, max_h] = max_halfplane_internal(drop_net, drop_red, drop_blue,
                                                     [](halfspace2_t const& h) { (void)h; return true;},
                                                     [&](double m, double b) { return f(m, m_Total, b, b_Total);}
            );

            //assert(util::aeq(evaluate_range(h, drop_red, drop_blue, f), max_h));
            if (max_h >= max_discrepancy) {
                max_discrepancy = max_h;
                max_halfspace = HalfSpace<3>(lift_half_space(h, pivot));
                //assert(util::aeq(evaluate_range(max_halfspace, red, blue, f), max_h));
            }

        }
        return std::make_tuple(halfspace3_t(max_halfspace), max_discrepancy);
    }


    std::tuple<halfspace3_t, double> max_halfspace_labeled(
            const point3_list_t& point_net,
            const lpoint3_list_t& red,
            const lpoint3_list_t& blue,
            const discrepancy_func_t& f) {
        return max_halfspace_internal(point_net, red, blue, f);
    }

    std::tuple<halfspace3_t, double> max_halfspace(
            const point3_list_t& point_net,
            const wpoint3_list_t& red,
            const wpoint3_list_t& blue,
            const discrepancy_func_t& f) {
        return max_halfspace_internal(point_net, red, blue, f);
    }



    template <typename T>
    void remove_duplicates_alt(T& pts) {
        std::sort(pts.begin(), pts.end(), [](auto const& p1, auto const& p2){
            return p1[0] < p2[0];

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


    lpoint_list_t compress_pts(lpoint_list_t lpoints) {


        auto comp = [] (lpt2_t const& p1, lpt2_t const& p2) {
            return p1.get_label() < p2.get_label();
        };
        std::sort(lpoints.begin(), lpoints.end(), comp) ;

        lpoint_list_t output_pts;
        for (auto pt_it = lpoints.begin(); pt_it != lpoints.end(); ) {
            auto end_of_group = std::upper_bound(pt_it, lpoints.end(), *pt_it, comp);

            //Split by c
            auto first_part = std::partition(pt_it, end_of_group, [](pt2_t const& pt) {
                return 0 < pt[2];
            });
            auto hull_1 = graham_march(point_list_t(pt_it, first_part));

            //Have to flip these signs
            point_list_t flipped_pts;
            for (auto it = first_part; it != end_of_group; it++) {
                flipped_pts.emplace_back(it->flip_orientation());
            }
            auto hull_2 = graham_march(flipped_pts);

            for (auto &p: hull_1) {
                output_pts.emplace_back(LPoint(pt_it->get_label(), pt_it->get_weight(), p[0], p[1], p[2]));
            }
            for (auto &p: hull_2) {
                //Flip these points back.
                output_pts.emplace_back(LPoint(pt_it->get_label(), pt_it->get_weight(), -p[0], -p[1], -p[2]));
            }
            pt_it = end_of_group;
        }
        return output_pts;

    }


    std::tuple<halfspace3_t, double> max_halfspace_restricted(
            const pt3_t& pivot,
            const lpoint3_list_t& point_net,
            const lpoint3_list_t& red,
            const lpoint3_list_t& blue,
            bool compress,
            const filter_func3_t& filter,
            const discrepancy_func2_t& f) {

        auto drop = [&pivot](const auto& pt) {
            return drop_point(pivot, pt);
        };
        auto transform_drop_net = [&pivot, &drop](const lpoint3_list_t& pts, lpoint_list_t& out) {
            //If the point is a duplicate of the pivot then we remove it
            for (auto &pt : pts) {
                if (!pt.approx_eq(pivot)) {
                    out.push_back(drop(pt));
                }
            }
        };

        auto transform_drop_lpoints = [&pivot, &drop](const lpoint3_list_t& pts, lpoint_list_t& out) {
            //If the point is a duplicate of the pivot then we remove all the lpoints corresponding to this value
            //and add this to the discrepancy_func2
            double removed_weight = 0;
            std::unordered_set<size_t> duplicate_labels;
            for (auto &pt : pts) {
                if (pt.approx_eq(pivot) && duplicate_labels.find(pt.get_label()) != duplicate_labels.end()) {
                    duplicate_labels.insert(pt.get_label());
                    removed_weight += pt.get_weight();
                }
            }
            for (auto &pt : pts) {
                if (duplicate_labels.find(pt.get_label()) == duplicate_labels.end()) {
                    out.push_back(drop(pt));
                }
            }
            return removed_weight;
        };



        lpoint_list_t drop_net_labeled;
        lpoint_list_t drop_red;
        lpoint_list_t drop_blue;
        transform_drop_net(point_net, drop_net_labeled);
        double add_red = transform_drop_lpoints(red, drop_red);
        double add_blue = transform_drop_lpoints(blue, drop_blue);

        if (compress) {
            drop_net_labeled = compress_pts(drop_net_labeled);
            drop_red = compress_pts(drop_red);
            drop_blue = compress_pts(drop_blue);
        }

        point_list_t drop_net;
        for (auto& pt : drop_net_labeled) {
            drop_net.push_back(pt);
        }
        auto func = [filter, pivot] (halfspace2_t const& h) {
            HalfSpace<3> halfspace_tmp(lift_half_space(h, pivot));
            return filter(halfspace_tmp);
        };

        auto [h, hmx] = max_halfplane_internal(drop_net, drop_red, drop_blue, func, [&](double m, double b){
            return f(m + add_red, b + add_blue);
        });
        return std::make_tuple(halfspace3_t(lift_half_space(h, pivot)), hmx);
    }

    std::tuple<halfspace3_t, double> max_halfspace_restricted(
            const pt3_t& pivot,
            const point3_list_t& point_net,
            const wpoint3_list_t& red,
            const wpoint3_list_t& blue,
            const filter_func3_t& filter,
            const discrepancy_func2_t& f) {

        auto drop = [&pivot](const auto& pt) {
            return drop_point(pivot, pt);
        };
        point_list_t drop_net(point_net.size(), pt2_t());
        wpoint_list_t drop_red(red.size(), wpt2_t());
        wpoint_list_t drop_blue(blue.size(), wpt2_t());
        std::transform(point_net.begin(), point_net.end(), drop_net.begin(), drop);
        std::transform(red.begin(), red.end(), drop_red.begin(), drop);
        std::transform(blue.begin(), blue.end(), drop_blue.begin(), drop);

        auto func = [filter, pivot] (halfspace2_t const& h) {
            HalfSpace<3> halfspace_tmp(lift_half_space(h, pivot));
            return filter(halfspace_tmp);
        };

        auto [h, hmx] = max_halfplane_internal(drop_net, drop_red, drop_blue, func, f);
        return std::make_tuple(halfspace3_t(lift_half_space(h, pivot)), hmx);
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


//    using seg_t = Segment;
//    using seg_list_t = std::vector<Segment>;
//
//    class CutCell {
//        double red_w = 0;
//        double blue_w = 0;
//        double cell_sc = 1.0;
//        seg_list_t red_segments;
//        seg_list_t blue_segments;
//    public:
//
//        CutCell() {}
//
//        CutCell(const point_list_t& red, const point_list_t& blue) {
//            for (auto const& pt : red) {
//                red_segments.emplace_back(pt);
//            }
//            for (auto const& pt : blue) {
//                blue_segments.emplace_back(pt);
//            }
//        }
//
//        CutCell(seg_list_t red, double rw, seg_list_t blue, double bw, double csc) :
//            red_w(rw), blue_w(bw), cell_sc(csc),
//            red_segments(std::move(red)), blue_segments(std::move(blue)) {}
//
//        double red_weight() const { return red_w; }
//        double blue_weight() const { return red_w; }
//        double cell_scale() const { return cell_sc; }
//
//        template <typename RNG>
//        std::tuple<CutCell, CutCell> split(RNG& rng) const {
//            // Uniformly
//
//            Segment seg;
//            std::uniform_real_distribution<double> dist(0, 1);
//            if ((blue_segments.empty() || dist(rng) < .5) && !red_segments.empty()) {
//                std::uniform_int_distribution<size_t> distribution(0, red_segments.size() - 1);
//                seg = red_segments[distribution(rng)];
//            } else if (!blue_segments.empty()){
//                std::uniform_int_distribution<size_t> distribution(0, blue_segments.size() - 1);
//                seg = blue_segments[distribution(rng)];
//            } else {
//                //Return two empty cells.
//                return {CutCell(), CutCell()};
//            }
//            auto split_segments = [&](seg_list_t const& crossing_segments) {
//
//                double below = 0;
//                seg_list_t upper_crossing_seg;
//                seg_list_t lower_crossing_seg;
//
//                for (auto& tseg : crossing_segments) {
//                    if (tseg.approx_eq(seg)) {
//                        below += cell_sc;
//                        continue;
//                    } else if (tseg.crossed(seg)) {
//                        auto [upper, lower] = tseg.split(seg);
//                        lower_crossing_seg.emplace_back(lower);
//                        upper_crossing_seg.emplace_back(upper);
//                    } else if (tseg.gte(seg)) {
//                        upper_crossing_seg.push_back(tseg);
//                    } else {
//                        lower_crossing_seg.push_back(tseg);
//                        below += cell_sc;
//                    }
//                }
//                return make_tuple(upper_crossing_seg, lower_crossing_seg, below);
//            };
//
//            auto [u_red_crossing, l_red_crossing, red_below_w] = split_segments(red_segments);
//            auto [u_blue_crossing, l_blue_crossing, blue_below_w] = split_segments(blue_segments);
//
//            //Subsample the lower and upper lines.
////            std::shuffle(u_red_crossing.begin(), u_red_crossing.end(), rng);
////            std::shuffle(l_red_crossing.begin(), l_red_crossing.end(), rng);
////            std::shuffle(u_blue_crossing.begin(), u_blue_crossing.end(), rng);
////            std::shuffle(l_blue_crossing.begin(), l_blue_crossing.end(), rng);
////            size_t new_cell_size = u_red_crossing.size() + l_red_crossing.size() +
////                    u_blue_crossing.size() + l_blue_crossing.size();
////
////            //double u_prop = (u_red_crossing.size() + u_blue_crossing.size()) / static_cast<double>(new_cell_size);
////            //double l_prop = 1 - u_prop;
////
////            double rescale = static_cast<double>(red_segments.size() + blue_segments.size()) / new_cell_size;
////
////            auto resize = [rescale](seg_list_t& pts) {
////                size_t new_size = lround(rescale * pts.size());
////                if (new_size >= pts.size()) {
////                    return ;
////                } else {
////                    pts.resize(new_size);
////                }
////            };
////
////            size_t ursz = u_red_crossing.size(), ubsz = u_blue_crossing.size(),
////                lrsz = l_red_crossing.size(), lbsz = l_blue_crossing.size();
////
////            resize(u_red_crossing);
////            resize(l_red_crossing);
////            resize(u_blue_crossing);
////            resize(l_blue_crossing);
////
////            double upper_scale = (ursz + ubsz) / static_cast<double>(u_red_crossing.size() + u_blue_crossing.size());
////            double lower_scale = (lrsz + lbsz) / static_cast<double>(l_red_crossing.size() + l_blue_crossing.size());
//            //Need the red lines and the blue lines to represent the original sets.
////            CutCell c1(u_red_crossing, red_below_w + red_w,
////                        u_blue_crossing, blue_below_w + blue_w,
////                        cell_sc * upper_scale);
////            CutCell c2(l_red_crossing, red_w, l_blue_crossing, blue_w, cell_sc * lower_scale);
//            CutCell c1(u_red_crossing, red_w + red_below_w,
//                       u_blue_crossing, blue_w + blue_below_w,
//                       1.0);
//            CutCell c2(l_red_crossing, red_w, l_blue_crossing, blue_w, 1.0);
//            return {c1, c2};
//        }
//
//        template<typename RNG>
//        std::tuple<Point<2>, double, double> choose_line(RNG& rng) const {
//            Segment seg;
//            std::uniform_real_distribution<double> dist(0, 1);
//            if ((blue_segments.empty() || dist(rng) < .5) && !red_segments.empty()) {
//                std::uniform_int_distribution<size_t> distribution(0, red_segments.size() - 1);
//                seg = red_segments[distribution(rng)];
//            } else if (!blue_segments.empty()){
//                std::uniform_int_distribution<size_t> distribution(0, blue_segments.size() - 1);
//                seg = blue_segments[distribution(rng)];
//            } else {
//                //Return the line at infinity
//                return {Point<2>(0.0, 0.0, 1.0), 0.0, 0.0};
//            }
//            //Choose an endpoint
//            auto line = seg.get_e1().orient_down();
//
//            double red_l_weight = red_w;
//            for (auto& pt : red_segments) {
//                if (line.above_closed(pt)) red_l_weight += cell_sc;
//            }
//
//            double blue_l_weight = blue_w;
//            for (auto& pt : blue_segments) {
//                if (line.above_closed(pt)) blue_l_weight += cell_sc;
//            }
//            return {line, red_l_weight, blue_l_weight};
//        }
//
//        double get_total_weight() const {
//            return cell_sc * (red_segments.size() + blue_segments.size());
//        }
//
//        bool operator<(const CutCell& s2) const {
//            return get_total_weight() < s2.get_total_weight();
//        }
//
//        seg_list_t get_red_set() const {
//            return red_segments;
//        }
//
//        seg_list_t get_blue_set() const {
//            return blue_segments;
//        }
//    };
//
//    std::tuple<halfspace2_t, double> max_halfplane_fast(size_t plane_count,
//                                                        const point_list_t &red,
//                                                        const point_list_t &blue,
//                                                        const discrepancy_func_t &f) {
//
//
//        auto red_total = static_cast<double>(red.size());
//        auto blue_total = static_cast<double>(blue.size());
//        std::random_device rd;
//        std::minstd_rand gen(rd());
//        std::priority_queue<CutCell> cells;
//
//        cells.emplace(red, blue);
//
//        while (cells.size() < plane_count) {
//
//            auto& curr_node = cells.top();
//
//            if (curr_node.get_total_weight() < 2) {
//                break;
//            }
//            auto [upper_cell, lower_cell] = curr_node.split(gen);
//            cells.pop();
//            cells.push(upper_cell);
//            cells.push(lower_cell);
//        }
//
//
//        //Now just choose a single vertex from each cell and compute its discrepancy.
//        double max_disc = 0;
//        HalfSpace<2> max_halfspace;
//        while (!cells.empty()) {
//            auto& curr_node = cells.top();
//            auto [line, red_w, blue_w] = curr_node.choose_line(gen);
//            double scan_val = f(red_w, red_total, blue_w, blue_total);
//            if (scan_val > max_disc) {
//                max_halfspace = HalfSpace<2>(line);
//                max_disc = scan_val;
//            }
//            cells.pop();
//        }
//        return {max_halfspace, max_disc};
//    }
}
