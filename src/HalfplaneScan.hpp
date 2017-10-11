//
// Created by mmath on 9/19/17.
//

#ifndef PYSCAN_HALFPLANE_HPP
#define PYSCAN_HALFPLANE_HPP

#include "Point.hpp"

namespace pyscan {

    template <int dim=2>
    class Halfspace {
        /*
         * Pretend the halfspace is of form 1 = a_1 x_1 + a_2 x_2 + ... + a_k x_k
         */
        double coords[dim] = {0};
        double value;

        template <typename F>
        Halfspace(int d, F val) {
            coords[dim - d] = val;
        }

        template <typename F, typename ...Coords>
        Halfspace(int d, F val, Coords ...rest) : Halfspace(d - 1, rest...) {

            coords[dim - d] = val;
        }

    public:
        Halfspace() : value(0) {}

        template <typename ...Coords>
        Halfspace(double val, Coords... args) : Halfspace(dim, args...) {
            value = val;
            static_assert(sizeof...(Coords) == dim, "coords must be same length as dim");
        }

        //template <typename ...Ps>
        //Halfspace(Halfspace<2> h1, Ps... args) {
            /*
             * Takes a halfspace in dimension
             */

        //}

        bool contains(Point<dim> const& pt) const {
            return dot<double, dim>(pt, coords) >= 1;
        }

        template <int ix>
        double get() const {
            static_assert(ix < dim && 0 <= ix, "ix has to be in [0, dim)");
            return coords[ix];
        }

        double fValue() const {
            return value;
        }

        template <int ix, int d> friend double get(Halfspace<d> const& h);
    };

    template <int ix, int d>
    double get(Halfspace<d> const& h) {
        return h.coords[ix];
    }

    Halfspace<3> liftHalfspace(Halfspace<2> const& h2, Point<3> const& p3);

    Point<2> dropPoint(Point<3> const& fixed_point, Point<3> const& p3);

    double angle(Point<> const & p1, Point<> const& p2);

    std::tuple<Halfspace<>, Point<>, Point<>>
    maxHalfplaneLin(point_it p_net_b, point_it p_net_e, point_it p_samp_b, point_it p_samp_e);

    std::tuple<Halfspace<>, Point<>, Point<>>
    maxHalfplaneStat(point_it p_net_b, point_it p_net_e, point_it p_samp_b, point_it p_samp_e, double rho);

    std::tuple<Halfspace<>, Point<>, Point<>>
    maxHalfplaneGamma(point_it p_net_b, point_it p_net_e, point_it p_samp_b, point_it p_samp_e, double rho);

}
#endif //PYSCAN_HALFPLANE_HPP
