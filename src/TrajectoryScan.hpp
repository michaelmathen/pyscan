//
// Created by mmath on 9/4/18.
//


#ifndef PYSCAN_TRAJECTORYSCAN_HPP
#define PYSCAN_TRAJECTORYSCAN_HPP

#include "RectangleScan.hpp"
#include "DiskScan.hpp"
#include "Point.hpp"

namespace pyscan {

    using point_it = point_list::iterator;

    struct traj_set {
        point_list traj_pts;
        std::vector<size_t> offsets;
        traj_set(point_list const& tp, std::vector<size_t> const& off) :
            traj_pts(tp),
            offsets(off){}
    };

    struct wtraj_set : public traj_set {
        std::vector<double> weights;
        wtraj_set(point_list const& tp,
                std::vector<size_t> const& off,
                std::vector<double> const& weights) : traj_set(tp, off), weights(weights) {}
    };




    //////////////////////////////////////
    //Full Scanning code//////////////////
    //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Disk scanning Trajectory code//////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::unordered_map<long, std::vector<Point<>>> grid_traj(point_it traj_b, point_it traj_e, double chord_l);

    std::unordered_map<long, std::vector<Point<>>> approximate_traj_cells(point_it traj_b,
            point_it traj_e,
            double chord_l,
            double eps);

    
    std::tuple<Disk, double> traj_disk_scan(traj_set &net,
                                            wtraj_set &sampleM,
                                            wtraj_set &sampleB,
                                            double alpha,
                                            double min_r,
                                            std::function<double(double, double)> const &scan);


    std::vector<Point<>> core_set_3d_traj(point_list traj, double alpha);

}

#endif //PYSCAN_TRAJECTORYSCAN_HPP
