//
// Created by mmath on 9/4/18.
//
#include <tuple>
#include <unordered_map>

#include "FunctionApprox.hpp"
#include "TrajectoryScan.hpp"
#include "TrajectoryCoreSet.hpp"
#include "DiskScan.hpp"

namespace pyscan {


    std::tuple<Disk, double> max_disk_traj_grid(trajectory_set_t const& net,
                                                 wtrajectory_set_t const& sampleM,
                                                 wtrajectory_set_t const& sampleB,
                                                 double alpha,
                                                 double max_r,
                                                 const discrepancy_func_t &scan) {

        double curr_r = alpha;
        double curr_disk_val = 0;
        Disk curr_disk;
        while (curr_r < max_r) {
            double chord_l = sqrt(4 * alpha * curr_r - 2 * alpha * alpha);
            // Go through each trajectory and approximate it with a coreset of points

            // Go through the set of net points.
            size_t label = 0;
            lpoint_list_t net_points;
            for (auto const &traj_set : net) {
                approx_traj_labels(traj_set.begin(), traj_set.end(), chord_l, alpha, label, 0.0, net_points);
                label++;
            }
            // go through the set of measured points
            lpoint_list_t sampleM_points;
            label = 0;
            for (auto const &b : sampleM) {
                approx_traj_labels(b.begin(), b.end(), chord_l, alpha, label, b.get_weight(), sampleM_points);
                label++;//increment label
            }
            lpoint_list_t sampleB_points;
            label = 0;
            for (auto const &b : sampleB) {
                approx_traj_labels(b.begin(), b.end(), chord_l, alpha, label, b.get_weight(), sampleB_points);
                label++;//increment label
            }
            // Scan the resulting set of points using standard labeled disk scanning function.

            // return the max disk.
            auto [disk, disk_val] = max_disk_scale_labeled(net_points, sampleM_points, sampleB_points, alpha, curr_r, scan);
            if (curr_disk_val < disk_val) {
               curr_disk = disk;
               curr_disk_val = disk_val;
            }
            curr_r *= 2;
        }
        return {curr_disk, curr_disk_val};
    }




}