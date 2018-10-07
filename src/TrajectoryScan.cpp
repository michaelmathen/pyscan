//
// Created by mmath on 9/4/18.
//
#include <tuple>
#include <unordered_map>

#include "FunctionApprox.hpp"
#include "TrajectoryScan.hpp"
#include "TrajectoryCoreSet.hpp"


namespace pyscan {





//    void core_set_3d_traj_internal(point_it traj_b, point_it traj_e, size_t label, double weight, double alpha,
//                          lpoint_list_t& output) {
//
//        auto approx = eps_core_set3(alpha, [&](Vec3 const& direction) {
//            auto pt_max = std::max_element(traj_b, traj_e, [&] (Point<> const& p1, Point<> const& p2){
//                double z_1 = p1(0) * p1(0) + p1(1) * p1(1);
//                double z_2 = p2(0) * p2(0) + p2(1) * p2(1);
//
//                double mag_1 = getX(p1) * direction[0] + getY(p1) * direction[1] + z_1 * direction[2];
//                double mag_2 = getX(p2) * direction[0] + getY(p2) * direction[1] + z_2 * direction[2];
//                return mag_1 < mag_2;
//            });
//            return Vec3{getX(*pt_max), getY(*pt_max), getX(*pt_max) * getX(*pt_max) + getY(*pt_max) * getY(*pt_max)};
//        });
//
//        for (auto& pt : approx) {
//            output.emplace_back(label, weight, pt[0], pt[1], 1.0);
//        }
//    }
    std::tuple<Disk, double> traj_disk_scan(trajectory_set_t const& net,
                                            wtrajectory_set_t const& sampleM,
                                            wtrajectory_set_t const& sampleB,
                                            double alpha,
                                            double min_r,
                                            std::function<double(double, double)> const &scan) {


        double chord_l = sqrt(4 * alpha * min_r - 2 * alpha * alpha);
        // Go through each trajectory and approximate it with a coreset of points

        // Go through the set of net points.
        point_list_t net_points;
        for(auto const& traj_set : net) {
            approx_traj(traj_set.begin(), traj_set.end(), chord_l, alpha, net_points);
        }
        // go through the set of measured points
        lpoint_list_t sampleM_points;
        size_t label = 0;
        for(auto const& b : sampleM) {
            approx_traj_labels(b.begin(), b.end(), chord_l, alpha, label, b.weight, sampleM_points);
            label++;//increment label
        }
        lpoint_list_t sampleB_points;
        label = 0;
        for(auto const& b : sampleB) {
            approx_traj_labels(b.begin(), b.end(), chord_l, alpha, label, b.weight, sampleB_points);
            label++;//increment label
        }
        // Scan the resulting set of points using standard labeled disk scanning function.

        // return the max disk.
        auto grid_r = static_cast<uint32_t>(1.0 / min_r);
        return disk_scan_scale(net_points, sampleM_points, sampleB_points, grid_r, scan);
    }

}