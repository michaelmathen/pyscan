//
// Created by mmath on 5/28/17.
//

#ifndef PYSCAN_POINT_HPP
#define PYSCAN_POINT_HPP
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <array>
#include <cassert>
#include <sstream>
#include <functional>

#include "Utilities.hpp"

//#include <boost/math/special_functions/next.hpp>

namespace pyscan {


    // Coordiates in homgenous space.
    template<int dim = 2>
    class Point {
    public:


        template <class ...Coords, std::enable_if_t<(sizeof...(Coords) == dim + 1)>* = nullptr>
        Point(Coords... rest) {

            static_assert(sizeof...(rest) == dim + 1, "Wrong number of arguments for the point type.");
            std::initializer_list<double> il({rest...});
            std::copy(il.begin(), il.end(), coords.begin());
        }

       /* Point<dim>& operator=(Point<dim> const& pt) {
            std::copy(coords.begin(), coords.end(), pt.coords.begin());
            return *this;
        }*/

        Point() {
            coords.fill(0.0);
        }

        friend std::ostream &operator<<(std::ostream &os, Point const &pt) {
            os << "Point(";
            for (auto &el: pt.coords) {
                os << el << ", ";
            }
            os << ")";
            return os;
        }

        double operator[](size_t i) const {
            assert(i < dim + 1);
            return coords[i];
        }

        double operator()(size_t i) const {
            assert(i < dim);
            return coords[i] / coords[dim];
        }

        std::string str() const {
            std::stringstream ss;
            ss << *this;
            return ss.str();
        }

        bool operator==(const Point<dim> &pt) const {
            for (size_t i = 0; i < dim + 1; ++i) {
                if (coords[i] != pt.coords[i])
                    return false;
            }
            return true;
        }

        inline double dot(const Point<dim> &other) const {
            double res = 0.0;
            for (size_t i = 0; i < dim + 1; ++i) {
                res += coords[i] * other.coords[i];
            }
            return res;
        }

        inline double square_dist(const Point<dim> &p) const {
            double res = 0.0;
            for (size_t i = 0; i < dim; ++i) {
                double tmp = (p[i] / p[dim] - coords[i] / coords[dim]);
                res += tmp * tmp;
            }
            return res;
        }

        inline bool approx_eq(Point<dim> const& p) const {
            return util::aeq(this->square_dist(p), 0);
        }

        inline bool above_closed(const Point<dim> &p) const {
            return util::alte(dot(p), 0.0);
        }

        inline bool below_closed(const Point<dim> &p) const {
            return util::alte(0.0, dot(p));
        }

        inline bool above(const Point<dim> &p) const {
            return util::alt(dot(p), 0.0);
        }

        inline bool below(const Point<dim> &p) const {
            return util::alt(0.0, dot(p));
        }

    private:
        std::array<double, dim + 1> coords;
    };

    template<int dim = 2>
    class WPoint : public Point<dim> {
    public:
        template<typename ...Coords>
        WPoint(double weight, Coords... rest)
                : Point<dim>(rest...), weight(weight) {}

        WPoint()
                : Point<dim>(), weight(0.0) {}

        inline double get_weight() const {
            return weight;
        }

    private:
        double weight;
    };


    template<int dim = 2>
    class LPoint : public WPoint<dim> {
    public:
        template<typename ...Coords>
        LPoint(size_t label, double weight, Coords... rest)
                : WPoint<dim>(weight, rest...), label(label) {}

        LPoint()
                : WPoint<dim>(), label(0) {}

        inline size_t get_label() const {
            return label;
        }

    private:
        size_t label;
    };


    Point<2> correct_orientation(const Point<2>& pivot, const Point<2>& p);

    Point<2> intersection(const Point<2> &p1, const Point<2> &p2);

    Point<3> cross_product(const Point<3> &p1, const Point<3> &p2);

    bool is_parallel(const Point<2> &l1, const Point<2> &l2);

    using pt2_t = Point<2>;
    using wpt2_t = WPoint<2>;
    using lpt2_t = LPoint<2>;
    using point_list_t = std::vector<pt2_t>;
    using point_it_t = point_list_t::iterator;
    using weight_list_t = std::vector<double>;
    using weight_it_t = weight_list_t::iterator;
    using pt3_t = Point<3>;
    using wpt3_t = WPoint<3>;
    using lpt3_t = LPoint<3>;
    using point3_list_t = std::vector<pt3_t>;
    using point3_it_t = point3_list_t::iterator;
    using label_list_t = std::vector<size_t>;
    using wpoint_list_t = std::vector<WPoint<2>>;
    using wpoint_it_t = wpoint_list_t::iterator;
    using lpoint_list_t = std::vector<LPoint<2>>;
    using lpoint_it_t = lpoint_list_t::iterator;
    using wpoint3_list_t = std::vector<WPoint<3>>;
    using lpoint3_list_t = std::vector<LPoint<3>>;
    using lpoint3_it_t = lpoint3_list_t::iterator;

    using discrepancy_func_t = std::function<double(double, double)>;

}


#endif //PYSCAN_POINT_HPP
