/*
 * Created by Michael Matheny on 6/3/19.
 * at the University of Utah
 * email: michaelmathen@gmail.com
 * website: https://mmath.dev/
 */
#include "Point.hpp"
#include "Disk.hpp"
#include "Gridding.hpp"
#include "SatScan.hpp"

namespace pyscan {

    static std::tuple<Disk, double> max_disk_sequence(
            const pt2_t& center,
            lpoint_list_t &measured,
            lpoint_list_t &baseline,
            discrepancy_func_t const& disc) {

        auto order_f = [&center] (const pt2_t& p1, const pt2_t& p2) {
            return center.dist(p1) < center.dist(p2);
        };
        std::sort(measured.begin(), measured.end(), order_f);
        std::sort(baseline.begin(), baseline.end(), order_f);
        double m_tot = computeTotal(measured);
        double b_tot = computeTotal(baseline);

        double max_disc = -std::numeric_limits<double>::infinity();
        Disk max_disk(center(0), center(1), 0.0);
        double m_curr_sum = 0;
        double b_curr_sum = 0;
        auto curr_m = measured.begin();
        auto curr_b = baseline.begin();
        double curr_dist = 0;
        using label_set = std::unordered_set<decltype(measured.front().get_label())>;
        label_set measured_label_set;
        label_set baseline_label_set;
        auto increment = [&] (const lpt2_t& pt, label_set& labels) {
            auto it = labels.find(pt.get_label());
            if (it == labels.end()) {
                labels.insert(pt.get_label());
                return pt.get_weight();
            } else {
                return 0.0;
            }
        };
        while (curr_m != measured.end() || curr_b != baseline.end()) {
            if (curr_b == baseline.end()) {
                m_curr_sum += increment(*curr_m, measured_label_set);
                curr_dist = center.dist(*curr_m);
                curr_m++;
            } else if (curr_m == measured.end()) {
                b_curr_sum += increment(*curr_b, baseline_label_set);
                curr_dist = center.dist(*curr_b);
                curr_b++;
            } else if (center.dist(*curr_m) < center.dist(*curr_b)) {
                m_curr_sum += increment(*curr_m, measured_label_set);
                curr_dist = center.dist(*curr_m);
                curr_m++;
            } else {
                b_curr_sum += increment(*curr_b, baseline_label_set);
                curr_dist = center.dist(*curr_b);
                curr_b++;
            }
            double curr_disc = disc(m_curr_sum, m_tot, b_curr_sum, b_tot);
            if (curr_disc > max_disc) {
                max_disc = curr_disc;
                max_disk = Disk(center(0), center(1), curr_dist);
            }
        }
        return std::make_tuple(max_disk, max_disc);
    }

    static std::tuple<Disk, double> max_disk_sequence(
            const pt2_t& center,
            wpoint_list_t &measured,
            wpoint_list_t &baseline,
            discrepancy_func_t const& disc) {

        auto order_f = [&center] (const pt2_t& p1, const pt2_t& p2) {
            return center.dist(p1) < center.dist(p2);
        };
        std::sort(measured.begin(), measured.end(), order_f);
        std::sort(baseline.begin(), baseline.end(), order_f);
        double m_tot = computeTotal(measured);
        double b_tot = computeTotal(baseline);

        double max_disc = -std::numeric_limits<double>::infinity();
        Disk max_disk(center(0), center(1), 0.0);
        double m_curr_sum = 0;
        double b_curr_sum = 0;
        auto curr_m = measured.begin();
        auto curr_b = baseline.begin();
        double curr_dist = 0;
        while (curr_m != measured.end() || curr_b != baseline.end()) {
            if (curr_b == baseline.end()) {
                m_curr_sum += curr_m->get_weight();
                curr_dist = center.dist(*curr_m);
                curr_m++;
            } else if (curr_m == measured.end()) {
                b_curr_sum += curr_b->get_weight();
                curr_dist = center.dist(*curr_b);
                curr_b++;
            } else if (center.dist(*curr_m) < center.dist(*curr_b)) {
                m_curr_sum += curr_m->get_weight();
                curr_dist = center.dist(*curr_m);
                curr_m++;
            } else {
                b_curr_sum += curr_b->get_weight();
                curr_dist = center.dist(*curr_b);
                curr_b++;
            }
            double curr_disc = disc(m_curr_sum, m_tot, b_curr_sum, b_tot);
            if (curr_disc > max_disc) {
                max_disc = curr_disc;
                max_disk = Disk(center(0), center(1), curr_dist);
            }
        }
       return std::make_tuple(max_disk, max_disc);
    }

    template <typename Pt>
    std::tuple<Disk, double> max_grid_disk_internal(
            std::vector<Pt> measured,
            std::vector<Pt> baseline,
            double grid_res,
            discrepancy_func_t const& disc,
            bbox_t const& full_bb) {
        Disk curr_max;
        double max_stat = -std::numeric_limits<float>::infinity();


        auto scanning = [&] (double x, double y){
            pt2_t center(x, y, 1.0);
            auto [d, mx_d] = max_disk_sequence(center, measured, baseline, disc);
            if (max_stat < mx_d) {
                curr_max = d;
                max_stat = mx_d;
            }

        };
        grid_scan(full_bb, grid_res, scanning);
        return std::make_tuple(curr_max, max_stat);
    }

    std::tuple<Disk, double> satscan_grid(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double disk_r,
            discrepancy_func_t const& func) {
        auto bb_op = bbox(measured, baseline);
        if (!bb_op.has_value()) {
            Disk curr_max;
            double max_stat = 0.0;
            return std::make_tuple(curr_max, max_stat);
        }
        auto full_bb = bb_op.value();
        auto [mnx, mny, mxx, mxy] = full_bb;
        auto edited_bb = std::make_tuple(mnx - disk_r, mny - disk_r, mxx + disk_r, mxy + disk_r);
        return max_grid_disk_internal(measured, baseline, grid_res, func, edited_bb);
    }

    std::tuple<Disk, double> satscan_grid_labeled(
            const lpoint_list_t &measured,
            const lpoint_list_t &baseline,
            double grid_res,
            double disk_r,
            discrepancy_func_t const& func) {
        auto bb_op = bbox(measured, baseline);
        if (!bb_op.has_value()) {
            Disk curr_max;
            double max_stat = 0.0;
            return std::make_tuple(curr_max, max_stat);
        }
        auto full_bb = bb_op.value();
        auto [mnx, mny, mxx, mxy] = full_bb;
        auto edited_bb = std::make_tuple(mnx - disk_r, mny - disk_r, mxx + disk_r, mxy + disk_r);
        return max_grid_disk_internal(measured, baseline, grid_res, func, edited_bb);
    }


}