/*
 * Created by Michael Matheny on 5/11/19.
 * at the University of Utah
 * email: michaelmathen@gmail.com
 * website: https://mmath.dev/
 */

#ifndef PYSCAN_GRIDDEDSCANNING_HPP
#define PYSCAN_GRIDDEDSCANNING_HPP

#include <tuple>
#include <vector>

#include "Point.hpp"
#include "Gridding.hpp"

/*
 *
if (net_chunk.size() >= dim) {
    for (auto pt1 = range.first; pt1 != range.second; ++pt1) {
        for (auto &pt2: net_chunk) {
            if (pt1->second.approx_eq(pt2)) continue;

            auto [local_max_disk, local_max_stat] = max_disk_restricted(pt1->second, pt2, net_chunk, red_chunk, blue_chunk,
                    min_res, 2 * min_res,
                    red_tot, blue_tot, f);
            if (local_max_stat > max_stat) {
                cur_max = local_max_disk;
                max_stat = local_max_stat;
            }

        }
    }
}
 */
/*
 * This is meant as a way to abstract out methods that uses a regular grid to prune out far away points. This
 * method was originally used in some of the disk scanning algorithms, but then became useful in a number of different
 * contexts. This is a nice method for converting
 */

template<typename Net_T, typename Point_T,  typename EvalF, typename Region_t>
std::tuple<Region_t, double> gridding_scan(
        const std::vector<Net_T> &point_net,
        const std::vector<Point_T> &red,
        const std::vector<Point_T> &blue,
        EvalF& eval,
        size_t boundary_width,
        double res) {
    Region_t curr_max;
    double max_stat = 0.0;
    if (point_net.empty()) {
        return std::make_tuple(curr_max, max_stat);
    }
    auto bb_op = bbox(point_net, red, blue);
    if (!bb_op.has_value()) {
        return std::make_tuple(curr_max, max_stat);
    }
    auto bb = bb_op.value();
    double red_tot = computeTotal(red);
    double blue_tot = computeTotal(blue);
    SparseGrid<Net_T> grid_net(bb, point_net, res);
    auto grid_r = grid_net.get_grid_size();
    SparseGrid<Point_T> grid_red(bb, red, res), grid_blue(bb, blue, res);


    for (auto center_cell = grid_net.begin(); center_cell != grid_net.end();) {
        std::vector<Net_T> net_chunk;
        std::vector<Point_T> red_chunk;
        std::vector<Point_T> blue_chunk;
        net_chunk.clear();
        red_chunk.clear();
        blue_chunk.clear();
        auto[i, j] = grid_net.get_cell(center_cell->second);
        size_t start_k = i < boundary_width ? 0 : i - boundary_width;
        size_t start_l = j < boundary_width ? 0 : j - boundary_width;
        size_t end_k = i + boundary_width < grid_r ? i + boundary_width : grid_r;
        size_t end_l = j + boundary_width < grid_r ? j + boundary_width : grid_r;
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

        auto[local_max_disk, local_max_stat] = eval(range, net_chunk, red_chunk, blue_chunk, red_tot, blue_tot);
        if (local_max_stat > max_stat) {
            curr_max = local_max_disk;
            max_stat = local_max_stat;
        }

        auto last = center_cell->first;
        do {
            ++center_cell;
        } while (center_cell != grid_net.end() && center_cell->first == last);
    }

    return std::make_tuple(curr_max, max_stat);
}

#endif //PYSCAN_GRIDDEDSCANNING_HPP
