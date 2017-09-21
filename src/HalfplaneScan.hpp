//
// Created by mmath on 9/19/17.
//

#ifndef PYSCAN_HALFPLANE_HPP
#define PYSCAN_HALFPLANE_HPP

#include "Point.hpp"

namespace pyscan {


    class Halfplane {
        double a, b;
        double value;
    public:
        Halfplane() : a(0), b(0), value(0) {}

        Halfplane(double a, double b, double val) : a(a), b(b), value(val) {}

        double getSlope() const {
            return a;
        }

        double getIntersect() const {
            return b;
        }

        double fValue() const {
            return value;
        }
    };

    double angle(Point<double> const & p1, Point<double> const& p2);

    Halfplane maxHalfplaneLin(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e);

    Halfplane maxHalfplaneStat(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e, double rho);

    Halfplane maxHalfplaneGamma(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e, double rho);
}
#endif //PYSCAN_HALFPLANE_HPP
