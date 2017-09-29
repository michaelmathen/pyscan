//
// Created by mmath on 9/25/17.
//

#ifndef PYSCAN_DISKSCAN_HPP_HPP
#define PYSCAN_DISKSCAN_HPP_HPP

#include <vector>

#include "Point.hpp"

namespace pyscan {
    class Disk {
        double value;
        double a;
        double b;
        double radius;
    public:
        Disk(double val, double a, double b, double r) : value(val), a(a), b(b), radius(r) {
        }

        Disk() : value(0), a(0), b(0), radius(0) {}

        template<typename W>
        bool contains(Point<W, 2> const &pt) const {
            return (get<0>(pt) - getA()) * (get<0>(pt) - getA()) +  (get<1>(pt) - getA()) * (get<1>(pt) - getA()) <= radius * radius;
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

        double fValue() const {
            return value;
        }
    };

    template <typename F>
    Disk diskScan(std::vector<Point<double, 2>>& net, std::vector<Point<double, 2>>& sampleM, std::vector<Point<double, 2>>& sampleB, F scan);

    Disk diskScanStatLabels(std::vector<LPoint<double, 2>>& net, std::vector<LPoint<double, 2>>& sampleM, std::vector<LPoint<double, 2>>& sampleB, double rho);



    Disk diskScanStat(std::vector<Point<double, 2>>& net, std::vector<Point<double, 2>>& sampleM, std::vector<Point<double, 2>>& sampleB, double rho);
}
#endif //PYSCAN_DISKSCAN_HPP_HPP
