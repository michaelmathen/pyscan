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

    template <int dim=2>
    class Halfspace {
        double coords[dim] = {0};
        double value;

        template <typename F>
        Halfspace(int d, F val) {
            coords[dim - d] = val;
            static_assert(std::is_same<double, F>::value, "coords must be double value");
        }

        template <typename F, typename ...Coords>
        Halfspace(int d, F val, Coords ...rest) : Halfspace(d - 1, rest...) {
            coords[dim - d] = val;
            static_assert(std::is_same<double, F>::value, "coords must be double value");
        }
    public:
        Halfspace() : value(0) {}

        template <typename ...Coords>
        Halfspace(double val, Coords... args) : Halfspace(dim, args...) {
            value = val;
            static_assert(sizeof...(Coords) == dim, "coords must be same length as dim");
        }

        template <int ix>
        double get() const {
            static_assert(ix < dim && 0 <= ix, "ix has to be in [0, dim)");
            return coords[ix];
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
