//
// Created by Michael on 10/5/2018.
//

#ifndef PYSCAN_TRAJECTORYCORESET_HPP
#define PYSCAN_TRAJECTORYCORESET_HPP

#include "Point.hpp"

namespace pyscan {


    point_list_t approx_traj_labels(point_list_t const& trajectory_pts, double weight, double chord_l, double eps);

    point3_list_t kernel3d(point3_list_t const& pts, double eps);
    point_list_t lifting_coreset(point_list_t const& pts, double eps);
}
#endif //PYSCAN_TRAJECTORYCORESET_HPP
