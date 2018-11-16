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
#include <unordered_map>

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

        friend Point<dim> operator*(Point<dim> pt, double val) {
            Point<dim> res;
            for (size_t i = 0; i < dim; ++i) {
                res.coords[i] = pt.coords[i] * val;
            }
            res.coords[dim] = pt.coords[dim];
            return res;
        }

        friend Point<dim> operator*(double val, Point<dim> pt) {
            return operator*(pt, val);
        }

        Point<dim> orient_down() const {
            // Always orient these down by the second to last coordinate.
            // In 2d this ensures that y is oriented down.
            // In 3d this ensures that z is oriented down.
            Point<dim> p_out;
            double orientation = std::copysign(1.0, coords[dim - 1]);
            for (size_t i = 0; i < dim + 1; i++) {
                p_out.coords[i] = -orientation * coords[i];
            }
            return p_out;
        }

        inline double& operator[](size_t i) {
            //assert(i < dim + 1);
            return coords[i];
        }

        inline double get_coord(size_t i) const {
            return coords[i];
        }

        inline double operator[](size_t i) const {
            //assert(i < dim + 1);
            return coords[i];
        }

        inline double operator()(size_t i) const {
            //assert(i < dim);
            return coords[i] / coords[dim];
        }

        virtual std::string str() const {
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

        inline double evaluate(const Point<dim> &other) const {
            double res = 0.0;
            for (size_t i = 0; i < dim + 1; ++i) {
                res += coords[i] * other.coords[i];
            }
            return res;
        }

        inline double norm_evaluate(const Point<dim> &other) const {
            /*
             * Insure that other is always a positive scaled value
             */
            double res = 0.0;
            for (size_t i = 0; i < dim + 1; ++i) {
                res += coords[i] * other.coords[i];
            }
            //double scale = std::copysign(1.0, other.coords[dim]);
            return res * 1.0;
        }

        inline double square_dist(const Point<dim> &p) const {
            double res = 0.0;
            for (size_t i = 0; i < dim; ++i) {
                double tmp = (p[i] / p[dim] - coords[i] / coords[dim]);
                res += tmp * tmp;
            }
            return res;
        }

        inline Point<dim> on_segment(Point<dim> const& p2, double alpha) {
            /*
             * Returns a point between this point and the second point that is alpha this + (1 - alpha) p2;
             */
            return alpha * (*this) + (1 - alpha) * p2;
        }

        inline double dist(const Point<dim> &p) const {
            return sqrt(square_dist(p));
        }


        Point<dim> operator-(const Point<dim>& other) const {
            Point<dim> res;
            for (size_t i = 0; i < dim; ++i) {
                res.coords[i] = coords[i] * other[dim] - other[i] * coords[dim];
            }
            res.coords[dim] = other[dim] * coords[dim];
            return res;
        }

        Point<dim> operator+(const Point<dim>& other) const {
            Point<dim> res;
            for (size_t i = 0; i < dim; ++i) {
                res.coords[i] = coords[i] * other[dim] + other[i] * coords[dim];
            }
            res.coords[dim] = other[dim] * coords[dim];
            return res;
        }


        inline double square_dist(const Point<dim>& begin, const Point<dim>& end) const {
            auto v = end - begin;
            auto w = *this - begin;
            double c1 = w.pdot(v);
            if (util::alte(c1, 0.0)) return square_dist(begin);
            double c2 = v.pdot(v);
            if (util::alte(c2, c1)) return square_dist(end);
            double b = c1 / c2;
            auto pb = begin + v * b;
            return square_dist(pb);
        }

        inline double pdot(const Point<dim>& other) const {
            /*
             * Takes the dot product between these points as if they are vectors
             */
            double res = 0.0;
            for (size_t i = 0; i < dim; ++i) {
                res += coords[i] * other.coords[i];
            }
            //If this is approximately 0 then the points exist at infinity and this operation doesn't make sense.
            assert(!util::aeq(coords[dim] * other[dim], 0.0));
            return res / (coords[dim] * other[dim]);
        }

        inline bool approx_eq(Point<dim> const& p) const {

            double res = 0.0;
            for (size_t i = 0; i < dim; ++i) {
                double tmp = (p[i] * coords[dim] - coords[i] * p[dim]);
                res += std::abs(tmp);
            }
            return util::aeq(res, 0.0);
        }

        inline bool above_closed(const Point<dim> &p) const {
            return util::alte(norm_evaluate(p), 0.0);
        }

        inline bool below_closed(const Point<dim> &p) const {
            return util::alte(0.0, norm_evaluate(p));
        }

        inline bool above(const Point<dim> &p) const {
            return util::alt(norm_evaluate(p), 0.0);
        }

        inline bool below(const Point<dim> &p) const {
            return util::alt(0.0, norm_evaluate(p));
        }


        inline bool crosses( const Point<dim>& p1, const Point<dim>& p2) const {
            /*
             * Checks if this line is crossed by this line segment between p1
             */
            return above(p1) != above(p2);
        }

        inline bool parallel_lte( const Point<dim>& l1) const {
            /*
             * Checks to see if this line is less than or equal to the line l1 assuming that l1 is parallel to
             * this.
             */
            double norm = 1.0;
            for (size_t i = 0; i < dim; ++i) {
                if (!util::aeq(l1[i], 0)) {
                    norm = coords[i] / l1[i];
                    break;
                }
            }
            return util::alte(l1[dim] * norm, coords[dim]);
        }

        inline bool parallel_lt( const Point<dim>& l1) const {
            /*
             * Checks to see if this line is less than or equal to the line l1 assuming that l1 is parallel to
             * this.
             */
            double norm = 1.0;
            for (size_t i = 0; i < dim; ++i) {
                if (!util::aeq(l1[i], 0)) {
                    norm = coords[i] / l1[i];
                    break;
                }
            }
            return util::alt(l1[dim] * norm, coords[dim]);
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

        virtual void set_weight(double w) {
            weight = w;
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

        virtual void set_label(size_t l) {
            label = l;
        }

    private:
        size_t label;
    };


    Point<2> correct_orientation(const Point<2>& pivot, const Point<2>& p);

    Point<2> intersection(const Point<2> &p1, const Point<2> &p2);

    std::tuple<double, double, double> normal(const Point<3> &p1, const Point<3> &p2, const Point<3> &p3);

    bool is_parallel(const Point<2> &l1, const Point<2> &l2);

    bool crosses_segment(const Point<2> &p1, const Point<2> &p2, const Point<2> &q1, const Point<2> &q2);

    //Return the two points on the line that are equidistance to some other point.
//    std::tuple<Point<2>, Point<2>> chord_pts(const Point<2> &line, const Point<2> &origin, double dist);

    using pt2_t = Point<2>;
    using wpt2_t = WPoint<2>;
    using lpt2_t = LPoint<2>;
    using point_list_t = std::vector<pt2_t>;
    using point_it_t = point_list_t::iterator;
    using cpoint_it_t = point_list_t::const_iterator;

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
    using cwpoint_it_t = wpoint_list_t::const_iterator;

    using lpoint_list_t = std::vector<LPoint<2>>;
    using lpoint_it_t = lpoint_list_t::iterator;
    using clpoint_it_t = lpoint_list_t::const_iterator;

    using wpoint3_list_t = std::vector<WPoint<3>>;
    using lpoint3_list_t = std::vector<LPoint<3>>;
    using lpoint3_it_t = lpoint3_list_t::iterator;

    using discrepancy_func_t = std::function<double(double, double, double, double)>;


    template <typename T>
    void remove_duplicates(T& pts) {
        std::sort(pts.begin(), pts.end(), [](pt2_t const& p1, pt2_t const& p2){
            return p1(0) < p2(1);
        });

        auto end_it = std::unique(pts.begin(), pts.end(), [] (pt2_t const& p1, pt2_t const& p2) {
            return p1.approx_eq(p2);
        });
        pts.erase(end_it, pts.end());
    }


}


#endif //PYSCAN_POINT_HPP
