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

    //for 2d it doesn't make sense to approximate a trajectory that has very few points in it.
    //since the compressed trajectory will probably have the same number of points. The algorithm also
    //breaks if you pass in two few points. This constant is to ensure this doesn't happen.

    inline Point<2> _initializer(Point<2> const& p1, Point<2> const& p2) {
        auto pt = intersection(p1, p2);
        double inv_norm = 1 / sqrt(pt[0] * pt[0] + pt[1] * pt[1]);
        return Point<2>(pt[0] * inv_norm, pt[1] * inv_norm, pt[2] * inv_norm).orient_down(1);
    }


    inline Point<3> _initializer(Point<3> const& p1, Point<3> const& p2, Point<3> const& p3) {
        auto [a, b, c] = normal(p1, p2, p3);
        double inv_norm = 1 / sqrt(a * a + b * b + c * c);
        double val = a * p3(0) + b * p3(1) + c * p3(2);
        auto plane = pt3_t{a * inv_norm, b * inv_norm, c * inv_norm, -val * inv_norm};
        //Orient down in the z-axis.
        return plane.orient_down(2);
    }


    template <int dim>
    class HalfSpace : public Range<dim> {
    protected:
        Point<dim> line;
    public:
        HalfSpace() : line() {}

        explicit HalfSpace(Point<dim> const& line) : line(line) {}

        template <class ...Coords, std::enable_if_t<(sizeof...(Coords) == dim)>* = nullptr>
        HalfSpace(Coords... pts) {
            line = _initializer(pts...);
        }

        bool contains(Point<dim> const& pt) const final {
            // ax + by + cz <= 0
            // a x / z + b * y / z + c >= 0/c
            //Want to ensure that points have positive c.
            return line.above_closed(pt);
        }

//        bool exact_contains(Point<dim> const& pt) const final {
//            // ax + by + cz <= 0
//            // a x / z + b * y / z + c >= 0/c
//            //Want to ensure that points have positive c.
//            return line.above_closed(pt);
//        }

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

        virtual std::string str() const {
            std::stringstream ss;
            ss << "HalfSpace(";
            ss << this->line;
            ss << ")";
            return ss.str();
        }
    };


    using filter_func2_t = std::function<bool(HalfSpace<2> const& pt)>;
    using filter_func3_t = std::function<bool(HalfSpace<3> const& pt)>;

    using halfspace2_t = HalfSpace<2>;
    using halfspace3_t = HalfSpace<3>;

    using discrepancy_func2_t = std::function<double(double, double)>;

    pt2_t drop_point(const pt3_t& fixed_point, const pt3_t& p);
    wpt2_t drop_point(const pt3_t& fixed_point, const wpt3_t& p);
    lpt2_t drop_point(const pt3_t& fixed_point, const lpt3_t& p);
    pt3_t lift_half_space(const halfspace2_t& h, const pt3_t& p);

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


//    std::tuple<halfspace2_t, double> max_halfplane_fast(size_t plane_count,
//            const point_list_t &red,
//            const point_list_t &blue,
//            const discrepancy_func_t &f);

    std::tuple<halfspace3_t, double> max_halfspace(
            const point3_list_t &point_net,
            const wpoint3_list_t &red,
            const wpoint3_list_t &blue,
            const discrepancy_func_t &f);

    std::tuple<halfspace3_t, double> max_halfspace_simple(
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

    std::tuple<halfspace3_t, double> max_halfspace_labeled_simple(
        const point3_list_t &point_net,
        const lpoint3_list_t &red,
        const lpoint3_list_t &blue,
        const discrepancy_func_t &f);

    std::tuple<halfspace3_t, double> max_halfspace_restricted(
        const pt3_t& pt,
        const lpoint3_list_t& point_net,
        const lpoint3_list_t& red,
        const lpoint3_list_t& blue,
        bool compress,
        const filter_func3_t& filter,
        const discrepancy_func2_t& f);

    std::tuple<halfspace3_t, double> max_halfspace_restricted(
        const pt3_t& pt,
        const point3_list_t& point_net,
        const wpoint3_list_t& red,
        const wpoint3_list_t& blue,
        const filter_func3_t& filter,
        const discrepancy_func2_t& f);

}

#endif

