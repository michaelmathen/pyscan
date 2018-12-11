//
// Created by mmath on 9/4/18.
//


#ifndef PYSCAN_TRAJECTORYSCAN_HPP
#define PYSCAN_TRAJECTORYSCAN_HPP

#include "RectangleScan.hpp"
#include "Disk.hpp"
#include "HalfSpaceScan.hpp"
#include "Point.hpp"
#include "Trajectory.hpp"

namespace pyscan {




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
                                                double max_r,
                                                const discrepancy_func_t &scan);

}
    


#endif //PYSCAN_TRAJECTORYSCAN_HPP
