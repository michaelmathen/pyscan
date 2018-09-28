//
// Created by mmath on 9/25/18.
//

#ifndef PYSCAN_HALFSPACESCAN_HPP
#define PYSCAN_HALFSPACESCAN_HPP

#include <tuple>
#include <functional>

#include "Point.hpp"


namespace pyscan {

    std::tuple<pt2_t, double> MaxHalfPlane(
            point_list_t& point_net,
            const wpoint_list_t& red,
            const wpoint_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<pt2_t, double> MaxHalfPlaneSimple(
            point_list_t& point_net,
            const wpoint_list_t& red,
            const wpoint_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<pt3_t, double> MaxHalfSpace(
            point3_list_t& point_net,
            const wpoint3_list_t& red,
            const wpoint3_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<pt2_t, double> MaxHalfPlaneLabeled(
            point_list_t& point_net,
            lpoint_list_t& red,
            lpoint_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<pt2_t, double> MaxHalfPlaneLabeledSimple(
            point_list_t& point_net,
            const lpoint_list_t& red,
            const lpoint_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<pt3_t, double>MaxHalfSpaceLabeled(
            point3_list_t& point_net,
            const lpoint3_list_t& red,
            const lpoint3_list_t& blue,
            const discrepancy_func_t& f);

}

#endif

