//
// Created by mmath on 9/4/18.
//


#ifndef PYSCAN_TRAJECTORYSCAN_HPP
#define PYSCAN_TRAJECTORYSCAN_HPP

#include "RectangleScan.hpp"
#include "DiskScan.hpp"
#include "Point.hpp"

namespace pyscan {


    using trajectory_t = point_list_t;

    class WTrajectory {
        double weight;
        trajectory_t trajectory_pts;

    public:
        WTrajectory(double w, trajectory_t pts) : weight(w), trajectory_pts(std::move(pts)) {}

        double get_weight() const {
            return weight;
        }

        cpoint_it_t begin() const {
            return trajectory_pts.begin();
        }

        cpoint_it_t end() const {
            return trajectory_pts.end();
        }

        point_it_t begin() {
            return trajectory_pts.begin();
        }

        point_it_t end() {
            return trajectory_pts.end();
        }
    };

    using wtrajectory_t = WTrajectory;
    using trajectory_set_t = std::vector<trajectory_t>;
    using wtrajectory_set_t = std::vector<wtrajectory_t>;


    //////////////////////////////////////
    //Full Scanning code//////////////////
    //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Disk scanning Trajectory code//////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::tuple<Disk, double> max_disk_traj_grid(trajectory_set_t const& net,
                                                wtrajectory_set_t const& sampleM,
                                                wtrajectory_set_t const& sampleB,
                                                double alpha,
                                                double min_r,
                                                const discrepancy_func_t &scan);
}
    


#endif //PYSCAN_TRAJECTORYSCAN_HPP
