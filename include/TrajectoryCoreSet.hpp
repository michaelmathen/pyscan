//
// Created by Michael on 10/5/2018.
//

#ifndef PYSCAN_TRAJECTORYCORESET_HPP
#define PYSCAN_TRAJECTORYCORESET_HPP

#include "Point.hpp"

namespace pyscan {

    void approx_traj(point_list_t::const_iterator traj_b, point_list_t::const_iterator traj_e,
                     double chord_l, double eps, point_list_t& output);

    void approx_traj_labels(point_list_t::const_iterator traj_b, point_list_t::const_iterator traj_e,
                            double chord_l, double eps, size_t label, double weight, lpoint_list_t& output);

    point_list_t approx_traj_labels(point_list_t const& trajectory_pts, double weight, double chord_l, double eps);

    point_list_t approx_traj_grid(point_list_t const& trajectory_pts, double chord_l, double eps);
    point3_list_t kernel3d(point3_list_t const& pts, double eps);
    point_list_t lifting_coreset(point_list_t const& pts, double eps);
}
#endif //PYSCAN_TRAJECTORYCORESET_HPP
