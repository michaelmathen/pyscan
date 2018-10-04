//
// Created by mmath on 9/25/18.
//

#ifndef PYSCAN_HALFSPACESCAN_HPP
#define PYSCAN_HALFSPACESCAN_HPP

#include <tuple>
#include <functional>

#include "Range.hpp"
#include "Point.hpp"


namespace pyscan {

    inline Point<2> _initializer(Point<2> const& p1, Point<2> const& p2) {
        return correct_orientation(p1, p2);
    }


    inline Point<3> _initializer(Point<3> const& p1, Point<3> const& p2, Point<3> const& p3) {
        auto normal = cross_product(p1, p2); //this is normalized in the last coordinate
        auto c = normal.dot(p3);
        return pt3_t{ normal(0), normal(1), normal(2), c};
    }


    template <int dim>
    class HalfSpace : public Range<dim> {
    private:
        Point<dim> line;
    public:
        HalfSpace() : line() {}

        explicit HalfSpace(Point<dim> const& line) : line(line) {}

        template <class ...Coords, std::enable_if_t<(sizeof...(Coords) == dim)>* = nullptr>
        HalfSpace(Coords... pts) {
            line = _initializer(pts...);
        }


        virtual bool contains(Point<dim> const& pt) const {
            return line.above_closed(pt);
        }

        auto get_coords() const -> decltype(line) {
            return line;
        }

        double operator[](size_t i) const {
            assert(i < dim + 1);
            return line[i];
        }

        double operator()(size_t i) const {
            assert(i < dim);
            return line(i);
        }


    };



    using halfspace2_t = HalfSpace<2>;
    using halfspace3_t = HalfSpace<3>;

    std::tuple<halfspace2_t, double> MaxHalfPlane(
            point_list_t& point_net,
            const wpoint_list_t& red,
            const wpoint_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<halfspace2_t, double> MaxHalfPlaneSimple(
            point_list_t& point_net,
            const wpoint_list_t& red,
            const wpoint_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<halfspace3_t, double> MaxHalfSpace(
            point3_list_t& point_net,
            const wpoint3_list_t& red,
            const wpoint3_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<halfspace2_t, double> MaxHalfPlaneLabeled(
            point_list_t& point_net,
            lpoint_list_t& red,
            lpoint_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<halfspace2_t, double> MaxHalfPlaneLabeledSimple(
            point_list_t& point_net,
            const lpoint_list_t& red,
            const lpoint_list_t& blue,
            const discrepancy_func_t& f);

    std::tuple<halfspace3_t, double>MaxHalfSpaceLabeled(
            point3_list_t& point_net,
            const lpoint3_list_t& red,
            const lpoint3_list_t& blue,
            const discrepancy_func_t& f);

}

#endif

