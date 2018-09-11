//
// Created by mmath on 9/4/18.
//


#ifndef PYSCAN_TRAJECTORYSCAN_HPP
#define PYSCAN_TRAJECTORYSCAN_HPP

#include "DiskScan.hpp"
#include "Point.hpp"

namespace pyscan {

    using point_it = point_list::iterator;

    struct traj_set {
        point_list traj_pts;
        std::vector<size_t> offsets;
        std::vector<double> weights;
    };

    struct wtraj_set : public traj_set {
        std::vector<double> weights;
    };


    /*
     * Compute the maximum disk region subject to the constraint a spatial closeness approximation.
     */
    std::tuple<Disk, double> traj_disk_scan(traj_set &net,
            wtraj_set &sampleM,
            wtraj_set &sampleB,
            double alpha,
            std::function<double(double, double)> const& scan);
}

#endif //PYSCAN_TRAJECTORYSCAN_HPP
