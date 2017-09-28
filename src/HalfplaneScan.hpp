//
// Created by mmath on 9/19/17.
//

#ifndef PYSCAN_HALFPLANE_HPP
#define PYSCAN_HALFPLANE_HPP

#include "Point.hpp"

namespace pyscan {


    template <typename W, int dim>
    double dot(pyscan::Point<W, dim> const& pt, double coords[dim]) {
        if (dim == 3) {
            return get<0>(pt) * coords[0] + get<1>(pt) * coords[1] + get<2>(pt) * coords[2];
        } else {
            return get<0>(pt) * coords[0] + get<1>(pt) * coords[1];
        }
    }



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

        template <typename W>
        bool contains(Point<W, dim> const& pt) const {
            return dot<W, dim>(pt, coords) >= 1;
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

    template <typename W>
    Halfspace<3> liftHalfspace(Halfspace<2> const& h2, Point<W, 3> const& p3) {
        return {h2.fValue(), get<0>(h2), get<1>(h2),
                (1 - get<0>(p3) * get<0>(h2) - get<1>(p3) * get<1>(h2)) / get<2>(p3)};
    }

    template <typename W>
    Point<W, 2> dropPoint(Point<W, 3> const& fixed_point, Point<W, 3> const& p3) {
        /*
         * Does an affine transformation from 3 to 2.
         */
        double scaling = get<2>(fixed_point) - get<2>(p3);
        return {p3.getRedWeight(), p3.getBlueWeight(),
                get<0>(p3) * get<2>(fixed_point) - get<0>(fixed_point) * get<2>(p3),
                get<1>(p3) * get<2>(fixed_point) - get<1>(fixed_point) * get<2>(p3)
        };
    }

    double angle(Point<double> const & p1, Point<double> const& p2);

    std::tuple<Halfspace<>, Point<double>, Point<double>>
    maxHalfplaneLin(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e);

    std::tuple<Halfspace<>, Point<double>, Point<double>>
    maxHalfplaneStat(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e, double rho);

    std::tuple<Halfspace<>, Point<double>, Point<double>>
    maxHalfplaneGamma(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e, double rho);

}
#endif //PYSCAN_HALFPLANE_HPP
