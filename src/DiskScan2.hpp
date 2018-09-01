//
// Created by mmath on 5/22/18.
//

#ifndef PYSCAN_DISKSCAN2_HPP
#define PYSCAN_DISKSCAN2_HPP

#include <functional>
#include "Point.hpp"

namespace pyscan {

    class Disk {
        double a;
        double b;
        double radius;
    public:
        Disk(double a_i, double b_i, double r) : a(a_i), b(b_i), radius(r) {
        }

        Disk() : a(0), b(0), radius(0) {}

        bool contains(Point<> const &pt) const {
            return (pt[0] / pt[2] - getA()) * (pt[0] / pt[2] - getA()) + (pt[1] / pt[2] - getB()) * (pt[1]/pt[2] - getB()) <= radius * radius;
        }

        void setA(double la) {
            a = la;
        }

        void setB(double lb) {
            b = lb;
        }

        double getA() const {
            return a;
        }

        double getB() const {
            return b;
        }

        double getR() const {
            return radius;
        }

        void setR(double rb) {
            radius = rb;
        }
    };

    std::tuple<Disk, double> max_disk(
            point_list& point_net,
            point_list& red,
            weight_list& red_w,
            point_list& blue,
            weight_list& blue_w,
            std::function<double(double, double)> const& f);

    std::tuple<Disk, double> max_disk_labeled(
            point_list& point_net,
            point_list& red,
            weight_list& red_w,
            label_list& red_labels,
            point_list& blue,
            weight_list& blue_w,
            label_list& blue_labels,
            std::function<double(double, double)> const& f);


//    std::tuple<Disk, double> max_disk_slow(point3_list& point_net,
//                       point3_list& red,
//                       weight_list& red_w,
//                       point3_list& blue,
//                       weight_list& blue_w,
//                       std::function<double(double, double)> const& f);
//
//    std::tuple<Disk, double> max_disk_slow_labeled(
//            point_list& point_net,
//            point_list& red,
//            weight_list& red_w,
//            label_list& red_labels,
//            point_list& blue,
//            weight_list& blue_w,
//            label_list& blue_labels,
//            std::function<double(double, double)> const& f);
}
#endif //PYSCAN_DISKSCAN2_HPP
