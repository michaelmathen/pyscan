//
// Created by mmath on 5/22/18.
//
#include <tuple>
#include <algorithm>
#include "DiskScan2.hpp"
#include "HalfSpaceScan.hpp"

namespace pyscan {


    lpt3_t llift_pt(lpt2_t const& pt) {
        double x = pt(0);
        double y = pt(1);
        return {pt.get_label(), pt.get_weight(), x, y, x * x + y * y, 1.0};
    }

    wpt3_t wlift_pt(wpt2_t const&  pt) {
        double x = pt(0);
        double y = pt(1);
        return {pt.get_weight(), x, y, x * x + y * y, 1.0};
    }

    pt3_t lift_pt(pt2_t const&  pt) {
        double x = pt(0);
        double y = pt(1);
        return {x, y, x * x + y * y, 1.0};
    }

    std::tuple<Disk, double> disk_scan2(
            point_list_t& point_net,
            wpoint_list_t& red,
            wpoint_list_t& blue,
            const discrepancy_func_t& f) {

        point3_list_t lifted_net(point_net.size(), pt3_t ());
        wpoint3_list_t lifted_red(red.size(), wpt3_t());
        wpoint3_list_t lifted_blue(blue.size(), wpt3_t());
        std::transform(point_net.begin(), point_net.end(), lifted_net.begin(), lift_pt);
        std::transform(red.begin(), red.end(), lifted_red.begin(), wlift_pt);
        std::transform(blue.begin(), blue.end(), lifted_blue.begin(), wlift_pt);
        auto mx_h = MaxHalfSpace(lifted_net, lifted_red, lifted_blue, f);
        auto h = std::get<0>(mx_h);
        double a = h[0], b = h[1], c = h[2], d = h[3];
        return std::make_tuple(Disk(-a / (2 * c), -b / (2 * c),
                                    sqrt((a * a + b * b - 4 * c * d) / (4 * c * c))),
                               std::get<1>(mx_h));
    };

    std::tuple<Disk, double> disk_scan2_labels(
            point_list_t& point_net,
            lpoint_list_t& red,
            lpoint_list_t& blue,
            const discrepancy_func_t& f) {

        point3_list_t lifted_net(point_net.size(), pt3_t());
        lpoint3_list_t lifted_red(red.size(), lpt3_t());
        lpoint3_list_t lifted_blue(blue.size(), lpt3_t ());
        std::transform(point_net.begin(), point_net.end(), lifted_net.begin(), lift_pt);
        std::transform(red.begin(), red.end(), lifted_red.begin(), llift_pt);
        std::transform(blue.begin(), blue.end(), lifted_blue.begin(), llift_pt);
        auto mx_h = MaxHalfSpaceLabeled(lifted_net, lifted_red,  lifted_blue, f);
        auto h = std::get<0>(mx_h);
        double a = h[0], b = h[1], c = h[2], d = h[3];
        return std::make_tuple(Disk(-a / (2 * c), -b / (2 * c),
                                    sqrt((a * a + b * b - 4 * c * d) / (4 * c * c))),
                               std::get<1>(mx_h));
    };
}