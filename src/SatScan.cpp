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

    std::tuple<Disk, double> max_grid_disk_internal(
            wpoint_list_t measured,
            wpoint_list_t baseline,
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


    std::tuple<Disk, double> satscan_points(
            wpoint_list_t measured,
            wpoint_list_t baseline,
            discrepancy_func_t const& func) {

        Disk curr_max;
        double max_stat = 0.0;

        for (auto &center : measured) {
            auto [d, mx_d] = max_disk_sequence(center, measured, baseline, func);
            if (max_stat < mx_d) {
                curr_max = d;
                max_stat = mx_d;
            }
        }
        for (auto& center: baseline) {
            auto [d, mx_d] = max_disk_sequence(center, measured, baseline, func);
            if (max_stat < mx_d) {
                curr_max = d;
                max_stat = mx_d;
            }
        }
        return std::make_tuple(curr_max, max_stat);
    }

    /*
     * This is meant as a comparison function for the kernel methods.
     */
//    std::tuple<Disk, double> max_grid_disk(
//            const wpoint_list_t &measured,
//            const wpoint_list_t &baseline,
//            double grid_res,
//            double radius_size) {
//
//
//        double red_tot = computeTotal(measured);
//        double blue_tot = computeTotal(baseline);
//        Disk curr_max;
//        double max_stat = 0.0;
//
//        auto bb_op = bbox(measured, baseline);
//        if (!bb_op.has_value()) {
//            return std::make_tuple(curr_max, max_stat);
//        }
//        auto full_bb = bb_op.value();
//
//        double x_offset_list[] = {0, radius_size, 0, radius_size};
//        double y_offset_list[] = {0, 0, radius_size, radius_size};
//        for (size_t o = 0; o < 4; o++) {
//
//            auto bb = std::make_tuple(
//                    x_offset_list[o] + std::get<0>(full_bb) - radius_size,
//                    std::get<1>(full_bb) + y_offset_list[o] - radius_size,
//                    std::get<2>(full_bb) + x_offset_list[o],
//                    std::get<3>(full_bb) + y_offset_list[o]);
//            //Define a coarse grid.
//            SparseGrid<wpt2_t> grid_red(bb, measured, radius_size * 3), grid_blue(bb, baseline, radius_size * 3);
//            std::vector<uint64_t> keys;
//            {
//                for (auto b = grid_red.begin(); b != grid_red.end(); b++) {
//                    keys.emplace_back(b->first);
//                }
//                for (auto b = grid_blue.begin(); b != grid_blue.end(); b++) {
//                    keys.emplace_back(b->first);
//                }
//                std::sort(keys.begin(), keys.end());
//                auto nend = std::unique(keys.begin(), keys.end());
//                keys.erase(nend, keys.end());
//            }
//
//            //Consider cells with points in them, so they must have $\eps n$ values
//            for (auto k : keys) {
//
//                auto[i, j] = grid_red.get_cell(k);
//
//                auto [lx, ly] = grid_red.get_lower_corner(k);
//                auto new_bb = std::make_tuple(lx + radius_size, ly + radius_size, lx + 2* radius_size, ly + 2 * radius_size);
//
//                wpoint_list_t measured_local;
//                wpoint_list_t baseline_local;
//                auto[r_it, r_e_it] = grid_red(i, j);
//
//                double subgrid_red_total = 0.0;
//                for (auto b = r_it; b != r_e_it; b++){
//                    subgrid_red_total += b->second.get_weight();
//                    measured_local.emplace_back(b->second);
//                }
//                auto[b_it, b_e_it] = grid_blue(i, j);
//                double subgrid_blue_total = 0.0;
//                for (auto b = b_it; b != b_e_it; b++){
//                    subgrid_blue_total += b->second.get_weight();
//                    baseline_local.emplace_back(b->second);
//                }
//                double dense_grid_res = grid_res * (red_tot + blue_tot) / (subgrid_blue_total + subgrid_red_total);
//
//                auto [max_d, max_v] = max_grid_disk_internal(
//                        measured_local,
//                        baseline_local,
//                        dense_grid_res,
//                        radius_size,
//                        disc, new_bb);
//                if (max_v > max_stat) {
//                    max_stat = max_v;
//                    curr_max = max_d;
//                }
//            }
//        }
//        return std::make_tuple(curr_max, max_stat);
//    }
}