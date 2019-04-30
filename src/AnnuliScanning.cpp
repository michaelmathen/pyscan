/*
 * Created by Michael Matheny on 4/25/19.
 * at the University of Utah
 * email: michaelmathen@gmail.com
 * website: https://mmath.dev/
 */

#include <vector>
#include <tuple>

#include "AnnuliScanning.hpp"
#include "SparseGrid.hpp"
#include "Utilities.hpp"


namespace pyscan {
    using disk_list_t = std::vector<Disk>;


    std::tuple<Disk, double> max_annuli(const point_list_t &pts,
                                        wpoint_list_t mpts,
                                        wpoint_list_t bpts,
                                        const std::vector<double> &radii,
                                        const discrepancy_func_t &func) {
        for (auto r_it = radii.begin(); r_it != radii.end(); ++r_it) {

            for (auto nit1 = pts.begin(); nit1 != pts.end() - 1; ++nit1) {
                for (auto nit2 = nit1 + 1; nit2 != pts.end(); ++nit2) {
                    // Check to make sure the points are close enough to be on the boundary of some disk of radii
                    // r
                    if (nit1->square_dist(*nit2) < 4 * (*r_it) * (*r_it)) {
                        //Now we can get the origin point and compute all the annuli
                        Disk test_disk(*nit1, *nit2, *r_it);
                        auto center = test_disk.getOrigin();

                        //Now we have the annuli
                        std::sort(mpts.begin(), mpts.end(), [&](const pt2_t& p1, const pt2_t& p2) {
                            return center.square_dist(p1) < center.square_dist(p2);
                        });

                        std::sort(bpts.begin(), bpts.end(), [&](const pt2_t& p1, const pt2_t& p2) {
                            return center.square_dist(p1) < center.square_dist(p2);
                        });
                        /*
                         * TODO write stuff here.
                         */
                    }
                }
            }
        }
    }

    std::tuple<Disk, double> max_annuli_restricted(
            pt2_t const& pt,
            const point_list_t &pts,
            wpoint_list_t mpts,
            wpoint_list_t bpts,
            const std::vector<double> &radii,
            double mtotal,
            double btotal,
            const discrepancy_func_t &func) {

        for (auto r_it = radii.begin(); r_it != radii.end(); ++r_it) {

            for (auto nit1 = pts.begin(); nit1 != pts.end() - 1; ++nit1) {
                // Check to make sure the points are close enough to be on the boundary of some disk of radii
                // r
                if (pt.square_dist(*nit1) < 4 * (*r_it) * (*r_it)) {
                    //Now we can get the origin point and compute all the annuli
                    Disk test_disk(*nit1, pt, *r_it);
                    auto center = test_disk.getOrigin();

                    //Now we have the annuli
                    std::sort(mpts.begin(), mpts.end(), [&](const pt2_t &p1, const pt2_t &p2) {
                        return center.square_dist(p1) < center.square_dist(p2);
                    });

                    std::sort(bpts.begin(), bpts.end(), [&](const pt2_t &p1, const pt2_t &p2) {
                        return center.square_dist(p1) < center.square_dist(p2);
                    });
                    /*
                     * TODO write stuff here.
                     */
                }
            }
        }

    }


    std::tuple<Disk, double> max_annuli_scale(
            const point_list_t &point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            const std::vector<double>& annuli_res,
            const discrepancy_func_t &f) {

        Disk cur_max;
        double max_stat = 0.0;
        if (point_net.empty() || annuli_res.empty()) {
            return std::make_tuple(Disk(), 0.0);
        }
        auto bb_op = bbox(point_net, red, blue);
        if (!bb_op.has_value()) {
            return std::make_tuple(cur_max, max_stat);
        }
        auto bb = bb_op.value();
        double red_tot = computeTotal(red);
        double blue_tot = computeTotal(blue);
        SparseGrid<pt2_t> grid_net(bb, point_net, annuli_res.front());
        auto grid_r = grid_net.get_grid_size();
        SparseGrid<wpt2_t> grid_red(bb, red, annuli_res.front()), grid_blue(bb, blue, annuli_res.front());


        for (auto center_cell = grid_net.begin(); center_cell != grid_net.end();) {
            std::vector<pt2_t> net_chunk;
            wpoint_list_t red_chunk;
            wpoint_list_t blue_chunk;
            net_chunk.clear();
            red_chunk.clear();
            blue_chunk.clear();
            size_t i, j;
            std::tie(i, j) = grid_net.get_cell(center_cell->second);
            size_t start_k = i < 2 ? 0 : i - 2;
            size_t start_l = j < 2 ? 0 : j - 2;
            size_t end_k = i + 2 < grid_r ? i + 2 : grid_r;
            size_t end_l = j + 2 < grid_r ? j + 2 : grid_r;
            auto range = grid_net(i, j);

            for (size_t k = start_k; k <= end_k; ++k) {
                for (size_t l = start_l; l <= end_l; ++l) {
                    auto net_range = grid_net(k, l);
                    for (auto it = net_range.first; it != net_range.second; ++it) {
                        net_chunk.emplace_back(it->second);
                    }

                    auto red_range = grid_red(k, l);
                    for (auto it = red_range.first; it != red_range.second; ++it)
                        red_chunk.emplace_back(it->second);


                    auto blue_range = grid_blue(k, l);
                    for (auto it = blue_range.first; it != blue_range.second; ++it)
                        blue_chunk.emplace_back(it->second);
                }
            }

            if (net_chunk.size() >= 3) {
                for (auto pt1 = range.first; pt1 != range.second; ++pt1) {
                    auto [local_max_disk, local_max_stat] =
                    max_annuli_restricted(pt1->second, net_chunk, red_chunk, blue_chunk,
                                        annuli_res,
                                        red_tot, blue_tot, f);
                    if (local_max_stat > max_stat) {
                        cur_max = local_max_disk;
                        max_stat = local_max_stat;
                    }
                }
            }


            auto last = center_cell->first;
            do {
                ++center_cell;
            } while (center_cell != grid_net.end() && center_cell->first == last);
        }

        return std::make_tuple(cur_max, max_stat);
    }

}