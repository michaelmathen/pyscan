//
// Created by mmath on 5/22/18.
//

#ifndef PYSCAN_DISKSCAN2_HPP
#define PYSCAN_DISKSCAN2_HPP

#include <functional>

#include "Point.hpp"
#include "Range.hpp"

namespace pyscan {


    class Disk : public Range<2> {

        pt2_t center;
        double radius;
    public:
        Disk(double a_i, double b_i, double r) : center(a_i, b_i, 1.0), radius(r) {
        }

        Disk(pt2_t const& p1, pt2_t const& p2, pt2_t const& p3);
        Disk() : center(0.0, 0.0, 1.0), radius(0) {}

        virtual bool contains(Point<> const &pt) const {
            return pt.square_dist(center) <= radius * radius;
        }

        double getA() const {
            return center(0);
        }

        double getB() const {
            return center(1);
        }

        double getR() const {
            return radius;
        }
    };

    std::tuple<Disk, double> disk_scan2(
            point_list_t& point_net,
            wpoint_list_t& red,
            wpoint_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<Disk, double> disk_scan2_labels(
            point_list_t& point_net,
            lpoint_list_t& red,
            lpoint_list_t& blue,
            const discrepancy_func_t& f);

}
#endif //PYSCAN_DISKSCAN2_HPP
