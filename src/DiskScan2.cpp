//
// Created by mmath on 5/22/18.
//
#include <tuple>
#include "DiskScan2.hpp"
#include "HalfplaneScan.hpp"

namespace pyscan {


    std::tuple<Pt3, double> max_disk(
            point_list& point_net,
            point_list& red,
            weight_list& red_w,
            point_list& blue,
            weight_list& blue_w,
            std::function<double(double, double)> const& f) {

        std::vector<Pt3> lifted_net(point_net.size(), Pt3());
        std::vector<Pt3> lifted_red(point_net.size(), Pt3());
        std::vector<Pt3> lifted_blue(point_net.size(), Pt3());
        auto lift = [] (Pt3 const& pt) {
            return Pt3(pt[0], pt[1], pt[0] * pt[0] + pt[1] * pt[1]);
        };
        std::transform(point_net.begin(), point_net.end(), lifted_net.begin(), lift);
        std::transform(red.begin(), red.end(), lifted_red.begin(), lift);
        std::transform(blue.begin(), blue.end(), lifted_blue.begin(), lift);
        auto mx_h = max_halfspace(lifted_net, lifted_red, red_w, lifted_blue, blue_w, f);
        auto h = std::get<0>(mx_h);
        double a = h[0], b = h[1], c = h[2], d = h[3];
        return { Pt3(-a / (2 * c), -b / (2 * c), (a * a + b * b - 4 * c * d) / (4 * c * c), 1), std::get<1>(mx_h)};
    };
}