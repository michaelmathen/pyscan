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
        auto pt = intersection(p1, p2).orient_down();
        double norm = sqrt(pt[0] * pt[0] + pt[1] * pt[1]);
        return Point<2>(pt[0] / norm, pt[1] / norm, pt[2] / norm);
    }


    inline Point<3> _initializer(Point<3> const& p1, Point<3> const& p2, Point<3> const& p3) {
        auto normal = cross_product(p1, p2); //this is normalized in the last coordinate
        auto c = normal.evaluate(p3);
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

        bool contains(Point<dim> const& pt) const final {
            return line.above_closed(pt);
        }

        bool intersects_segment(Point<dim> const& p1, Point<dim> const& p2) const final {
            return contains(p1) || contains(p2);
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

    std::tuple<halfspace2_t, double> max_halfplane(
            const point_list_t& point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            const discrepancy_func_t &f);

    std::tuple<halfspace2_t, double> max_halfplane_simple(
            const point_list_t &point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            const discrepancy_func_t &f);


    std::tuple<halfspace2_t, double> max_halfplane_fast(size_t plane_count,
            const point_list_t &red,
            const point_list_t &blue,
            const discrepancy_func_t &f);

    std::tuple<halfspace3_t, double> max_halfspace(
            const point3_list_t &point_net,
            const wpoint3_list_t &red,
            const wpoint3_list_t &blue,
            const discrepancy_func_t &f);

    std::tuple<halfspace2_t, double> max_halfplane_labeled(
            const point_list_t& point_net,
            const lpoint_list_t& red,
            const lpoint_list_t& blue,
            const discrepancy_func_t &f);

    std::tuple<halfspace2_t, double> max_halfplane_labeled_simple(
            const point_list_t &point_net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            const discrepancy_func_t &f);

    std::tuple<halfspace3_t, double> max_halfspace_labeled(
            const point3_list_t &point_net,
            const lpoint3_list_t &red,
            const lpoint3_list_t &blue,
            const discrepancy_func_t &f);

}

#endif

