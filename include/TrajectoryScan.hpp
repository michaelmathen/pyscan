//
// Created by mmath on 9/4/18.
//


#ifndef PYSCAN_TRAJECTORYSCAN_HPP
#define PYSCAN_TRAJECTORYSCAN_HPP

#include "RectangleScan.hpp"
#include "DiskScan.hpp"
#include "Point.hpp"

namespace pyscan {


    struct trajectory {
        point_list_t pts;

        virtual point_list_t::const_iterator begin() const {
            return pts.begin();
        }

        virtual point_list_t::const_iterator end() const {
            return pts.end();
        }
    };

    struct wtrajectory : public trajectory {
        double weight;
    };

    using trajectory_set_t = std::vector<trajectory>;
    using wtrajectory_set_t = std::vector<wtrajectory>;


    //////////////////////////////////////
    //Full Scanning code//////////////////
    //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Disk scanning Trajectory code//////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::tuple<Disk, double> traj_multilevel_disk_scan(trajectory_set_t const &net,
                                                       wtrajectory_set_t const &sampleM,
                                                       wtrajectory_set_t const &sampleB,
                                                       double alpha,
                                                       double min_r,
                                                       std::function<double(double, double)> const &scan);
}
    


#endif //PYSCAN_TRAJECTORYSCAN_HPP
