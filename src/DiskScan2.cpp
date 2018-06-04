//
// Created by mmath on 5/22/18.
//
#include <tuple>
#include "DiskScan2.hpp"
#include "HalfplaneScan.hpp"

namespace pyscan {


    std::tuple<Disk, double> max_disk(
            point_list& point_net,
            point_list& red,
            weight_list& red_w,
            point_list& blue,
            weight_list& blue_w,
            std::function<double(double, double)> const& f) {

        std::vector<Pt3> lifted_net(point_net.size(), Pt3());
        std::vector<Pt3> lifted_red(red.size(), Pt3());
        std::vector<Pt3> lifted_blue(blue.size(), Pt3());
        auto lift = [] (Pt2 const& pt) {
            double x = pt[0] / pt[2];
            double y = pt[1] / pt[2];
            return Pt3(x, y, x * x + y * y, 1.0);
        };
        std::transform(point_net.begin(), point_net.end(), lifted_net.begin(), lift);
        std::transform(red.begin(), red.end(), lifted_red.begin(), lift);
        std::transform(blue.begin(), blue.end(), lifted_blue.begin(), lift);
        auto mx_h = max_halfspace(lifted_net, lifted_red, red_w, lifted_blue, blue_w, f);
        auto h = std::get<0>(mx_h);
        double a = h[0], b = h[1], c = h[2], d = h[3];
        return std::make_tuple(Disk(-a / (2 * c), -b / (2 * c),
                                    sqrt((a * a + b * b - 4 * c * d) / (4 * c * c))),
                               std::get<1>(mx_h));
    };

    std::tuple<Disk, double> max_disk_labeled(
            point_list& point_net,
            point_list& red,
            weight_list& red_w,
            label_list& red_labels,
            point_list& blue,
            weight_list& blue_w,
            label_list& blue_labels,
            std::function<double(double, double)> const& f) {

        std::vector<Pt3> lifted_net(point_net.size(), Pt3());
        std::vector<Pt3> lifted_red(red.size(), Pt3());
        std::vector<Pt3> lifted_blue(blue.size(), Pt3());
        auto lift = [] (Pt2 const& pt) {
            double x = pt[0] / pt[2];
            double y = pt[1] / pt[2];
            return Pt3(x, y, x * x + y * y, 1.0);
        };
        std::transform(point_net.begin(), point_net.end(), lifted_net.begin(), lift);
        std::transform(red.begin(), red.end(), lifted_red.begin(), lift);
        std::transform(blue.begin(), blue.end(), lifted_blue.begin(), lift);
        auto mx_h = max_halfspace_labeled(lifted_net, lifted_red, red_w, red_labels, lifted_blue, blue_w, blue_labels, f);
        auto h = std::get<0>(mx_h);
        double a = h[0], b = h[1], c = h[2], d = h[3];
        return std::make_tuple(Disk(-a / (2 * c), -b / (2 * c),
                                    sqrt((a * a + b * b - 4 * c * d) / (4 * c * c))),
                               std::get<1>(mx_h));
    };
}