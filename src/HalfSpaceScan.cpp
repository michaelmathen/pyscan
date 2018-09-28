#include <algorithm>

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

inline static double under_line_labeled(const pt2_t& line, const lpoint_list_t& pts) {
    std::vector<bool> vis(pts.size(), false);
    double res = 0.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (line.above_closed(pts[i]) && !vis[pts[i].get_label()]) {
            res += pts[i].get_weight();
            vis[pts[i].get_label()] = true;
        }
    }
    return res;
}

inline static double under_line(const pt2_t& line, const wpoint_list_t& pts) {
    double res = 0.0;
    for (size_t i = 0; i < pts.size(); ++i)
        if (line.above_closed(pts[i]))
            res += pts[i].get_weight();
    return res;
}

std::tuple<pt2_t, double> MaxHalfPlane(
        point_list_t& point_net,
        const wpoint_list_t& red,
        const wpoint_list_t& blue,
        const discrepancy_func_t& f) {

    double red_total = sum_weight(red);
    double blue_total = sum_weight(blue);

    double max_discrepancy = -std::numeric_limits<double>::infinity();
    pt2_t max_line;

    assert(point_net.size() >= 2);
    for (size_t i = 0; i < point_net.size() - 1; ++i) {
        auto pivot = point_net[i];
        std::sort(point_net.begin() + i + 1, point_net.end(),
                  [&](const pt2_t& p1, const pt2_t& p2) { return calc_angle(pivot, p1) < calc_angle(pivot, p2); });

        auto l1 = correct_orientation(pivot, point_net[i + 1]);
        auto left_size = point_net.size() - i - 1;
        std::vector<double> red_delta(left_size, 0.0);
        std::vector<double> blue_delta(left_size, 0.0);
        std::vector<double> angles(left_size, 0.0);
        for (size_t j = 0; j < left_size; ++j)
            angles[j] = calc_angle(pivot, point_net[i + j + 1]);
        
        auto calc_delta = [&](const wpoint_list_t& pts, std::vector<double>& deltas) {
            double res = 0.0;
            for (size_t i = 0; i < pts.size(); ++i) {
                auto angle_it = std::upper_bound(angles.begin(), angles.end(), calc_angle(pivot, pts[i]));
                if (angle_it == angles.begin()) {
                    if (l1.above_closed(pts[i])) {
                        res += pts[i].get_weight();
                    }
                } else {
                    size_t ix = angle_it - angles.begin() - 1;
                    if (l1.above_closed(pts[i])) {
                        deltas[ix] -= pts[i].get_weight();
                        res += pts[i].get_weight();
                    } else {
                        deltas[ix] += pts[i].get_weight();
                    }
                }
            }
            return res;
        };

        double red_curr = calc_delta(red, red_delta);
        double blue_curr = calc_delta(blue, blue_delta);
        for (size_t j = 0; j < left_size; ++j) {
            double stat = f(red_curr / red_total, blue_curr / blue_total);
            if (max_discrepancy <= stat) {
                max_line = correct_orientation(pivot, point_net[i + j + 1]);
                max_discrepancy = stat;
            }
            red_curr += red_delta[j];
            blue_curr += blue_delta[j];
        }
    }

    return std::make_tuple(max_line, max_discrepancy);
}


std::tuple<pt2_t, double> MaxHalfPlaneSimple(
        point_list_t& point_net,
        const wpoint_list_t& red,
        const wpoint_list_t& blue,
        const discrepancy_func_t& f) {

    double max_discrepancy = -std::numeric_limits<double>::infinity();
    pt2_t max_line;
    
    double red_total = sum_weight(red);
    double blue_total = sum_weight(blue);
    assert(point_net.size() >= 2);
    for (size_t i = 0; i < point_net.size() - 1; ++i) {
        for (size_t j = i + 1; j < point_net.size(); ++j) {
            auto line = correct_orientation(point_net[i], point_net[j]);
            double stat = f(under_line(line, red) / red_total, under_line(line, blue) / blue_total);
            if (max_discrepancy <= stat) {
                max_line = line;
                max_discrepancy = stat;
            }
        }
    }
    return std::make_tuple(max_line, max_discrepancy);
}

std::tuple<pt3_t, double> MaxHalfSpace(
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
        auto pair_mx = MaxHalfPlane(drop_net, drop_red, drop_blue, f);
        if (std::get<1>(pair_mx) >= max_discrepancy) {
            max_discrepancy = std::get<1>(pair_mx);
            max_halfspace = lift_half_space(std::get<0>(pair_mx), pivot);
        }
    }
    return std::make_tuple(max_halfspace, max_discrepancy);
}

std::tuple<pt2_t, double> MaxHalfPlaneLabeled(
        point_list_t& point_net,
        lpoint_list_t& red,
        lpoint_list_t& blue,
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
            double stat = f(red_curr / red_total, blue_curr / blue_total);
            if (max_discrepancy <= stat) {
                max_line = correct_orientation(pivot, point_net[i + j + 1]);
                max_discrepancy = stat;
            }
            red_curr += red_delta[j];
            blue_curr += blue_delta[j];
        }
    }
    return std::make_tuple(max_line, max_discrepancy);
}


std::tuple<pt2_t, double> MaxHalfPlaneLabeledSimple(
        point_list_t& point_net,
        const lpoint_list_t& red,
        const lpoint_list_t& blue,
        const discrepancy_func_t& f) {

    pt2_t max_line;
    double max_discrepancy = -std::numeric_limits<double>::infinity();

    double red_total = sum_weight_unique(red);
    double blue_total = sum_weight_unique(blue);

    for (size_t i = 0; i < point_net.size() - 1; ++i) {
        for (size_t j = i + 1; j < point_net.size(); ++j) {
            auto line = correct_orientation(point_net[i], point_net[j]);
            double stat = f(under_line_labeled(line, red) / red_total, under_line_labeled(line, blue) / blue_total);
            if (max_discrepancy <= stat) {
                max_line = line;
                max_discrepancy = stat;
            }
        }
    }
    return std::make_tuple(max_line, max_discrepancy);
}

std::tuple<pt3_t, double> MaxHalfSpaceLabeled(
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
        auto pair_mx = MaxHalfPlaneLabeled(drop_net, drop_red, drop_blue, f);
        if (std::get<1>(pair_mx) >= max_discrepancy) {
            max_discrepancy = std::get<1>(pair_mx);
            max_halfspace = lift_half_space(std::get<0>(pair_mx), pivot);
        }
    }
    return std::make_tuple(max_halfspace, max_discrepancy);
}

}
