#include <algorithm>
#include <queue>
#include <random>

#include "Range.hpp"
#include "Segment.hpp"
#include "HalfSpaceScan.hpp"

namespace pyscan{




    inline static double calc_angle(const pt2_t& p1, const pt2_t& p2) {
        double y = p1[2] * p2[0] - p2[2] * p1[0];
        double x = p1[2] * p2[1] - p2[2] * p1[1];
        if (util::alte(0.0, y)) return atan2(y, x);
        else return M_PI + atan2(y, x);
    }



    inline static pt2_t drop_point(const pt3_t& fixed_point, const pt3_t& p) {
        return {p[0] * fixed_point[2] - fixed_point[0] * p[2],
                p[1] * fixed_point[2] - fixed_point[1] * p[2],
                p[3] * fixed_point[2] - fixed_point[3] * p[2]};
    }

    inline static wpt2_t drop_wpoint(const pt3_t& fixed_point, const wpt3_t& p) {
        return {p.get_weight(),
                p[0] * fixed_point[2] - fixed_point[0] * p[2],
                p[1] * fixed_point[2] - fixed_point[1] * p[2],
                p[3] * fixed_point[2] - fixed_point[3] * p[2]};
    }

    inline static lpt2_t drop_lpoint(const pt3_t& fixed_point, const lpt3_t& p) {
        return {p.get_label(),
                p.get_weight(),
                p[0] * fixed_point[2] - fixed_point[0] * p[2],
                p[1] * fixed_point[2] - fixed_point[1] * p[2],
                p[3] * fixed_point[2] - fixed_point[3] * p[2]};
    }

    inline static pt3_t lift_half_space(const pt2_t& h, const pt3_t& p) {
        return {h[0] * p[2],
                h[1] * p[2],
                -p[0] * h[0] - p[1] * h[1] - p[3] * h[2],
                h[2] * p[2]};
    }

    inline static double sum_weight(const wpoint_list_t& pts) {
        double res = 0.0;
        for (auto& x: pts) res += x.get_weight();
        return res;
    }

    inline static double sum_weight_unique(const lpoint_list_t& pts) {
        std::vector<bool> vis(pts.size(), false);
        double res = 0.0;
        for (auto& x: pts) {
            if (!vis[x.get_label()]) {
                res += x.get_weight();
                vis[x.get_label()] = true;
            }
        }
        return res;
    }


    std::tuple<halfspace2_t, double> max_halfplane(
            const point_list_t& point_net,
            const wpoint_list_t& red,
            const wpoint_list_t& blue,
            const discrepancy_func_t& f) {

        double red_total = sum_weight(red);
        double blue_total = sum_weight(blue);

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
                    halfspace2_t plane(pivot, pt);
                    auto angle_it = std::lower_bound(angles.begin(), angles.end(), -plane[0]);
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
            point3_list_t& point_net,
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

    std::tuple<halfspace2_t, double> max_halfplane_labeled(
            point_list_t point_net,
            lpoint_list_t red,
            lpoint_list_t blue,
            const discrepancy_func_t& f) {

        pt2_t max_line;
        double max_discrepancy = -std::numeric_limits<double>::infinity();

        double red_total = sum_weight_unique(red);
        double blue_total = sum_weight_unique(blue);

        for (size_t i = 0; i < point_net.size() - 1; ++i) {
            auto pivot = point_net[i];
            std::sort(point_net.begin() + i + 1, point_net.end(),
                      [&](const pt2_t& p1, const pt2_t& p2) { return calc_angle(pivot, p1) < calc_angle(pivot, p2); });
            auto l1 = correct_orientation(pivot, point_net[i + 1]);
            auto left_size = point_net.size() - i - 1;
            std::vector<double> red_delta(left_size, 0.0);
            std::vector<double> blue_delta(left_size, 0.0);

            auto calc_delta = [&](lpoint_list_t& pts, std::vector<double>& deltas) {
                std::vector<size_t> label_counts(pts.size(), 0);
                double res = 0.0;
                for (auto& x: pts) {
                    if (l1.above_closed(x)) {
                        if (label_counts[x.get_label()] == 0)
                            res += x.get_weight();
                        ++label_counts[x.get_label()];
                    }
                }

                std::sort(pts.begin(), pts.end(),
                        [&](const lpt2_t& x, const lpt2_t& y) { return calc_angle(pivot, x) < calc_angle(pivot, y); });

                size_t k = 0;
                size_t j = i + 1;
                while (k < pts.size() && calc_angle(pivot, pts[k]) < calc_angle(pivot, point_net[j])) ++k;
                for (++j; j < point_net.size(); ++j) {
                    double delta = 0.0;
                    for (; k < pts.size() && calc_angle(pivot, pts[k]) < calc_angle(pivot, point_net[j]); ++k) {
                        if (l1.above_closed(pts[k])) {
                            --label_counts[pts[k].get_label()];
                            if (label_counts[pts[k].get_label()] <= 0) {
                                delta -= pts[k].get_weight();
                            }
                        } else {
                            if (label_counts[pts[k].get_label()] == 0) {
                                delta += pts[k].get_weight();
                            }
                            ++label_counts[pts[k].get_label()];
                        }
                    }
                    deltas[j - i - 2] = delta;
                }

                return res;
            };

            double red_curr = calc_delta(red, red_delta);
            double blue_curr = calc_delta(blue, blue_delta);
            for (size_t j = 0; j < left_size; ++j) {
                double stat = f(red_curr, red_total, blue_curr, blue_total);
                if (max_discrepancy <= stat) {
                    max_line = correct_orientation(pivot, point_net[i + j + 1]);
                    max_discrepancy = stat;
                }
                red_curr += red_delta[j];
                blue_curr += blue_delta[j];
            }
        }
        return std::make_tuple(halfspace2_t(max_line), max_discrepancy);
    }


    std::tuple<halfspace3_t, double> max_halfspace_labeled(
            point3_list_t& point_net,
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


    std::tuple<halfspace2_t, double> max_halfplane_simple(
            point_list_t &point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            const discrepancy_func_t &f) {
        return max_range2<halfspace2_t, 2>(point_net, red, blue, f);
    }


    std::tuple<halfspace2_t, double> max_halfplane_labeled_simple(
            point_list_t &point_net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            const discrepancy_func_t &f) {
        return max_range2_labeled<halfspace2_t, 2>(point_net, red, blue, f);
    }

    std::tuple<halfspace3_t, double> max_halfspace_labeled_simple(
            point3_list_t &point_net,
            const lpoint3_list_t &red,
            const lpoint3_list_t &blue,
            const discrepancy_func_t &f) {
        return max_range3_labeled<halfspace3_t, 3>(point_net, red, blue, f);
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
