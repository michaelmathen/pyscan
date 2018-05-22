//
// Created by mmath on 5/28/17.
//

#ifndef PYSCAN_POINT_HPP
#define PYSCAN_POINT_HPP
#include <stdexcept>
#include <string>
#include <sstream>
#include <array>
#include <boost/math/special_functions/next.hpp>

namespace pyscan {

    using namespace boost::math;

#define MAX_FLOAT_DISTANCE 5

    inline bool aeq(double a, double b) {
        return fabs(float_distance(a, b)) < MAX_FLOAT_DISTANCE;
    }

    inline bool alt(double a, double b) {
        return a < b && !aeq(a, b);
    }

    inline bool alte(double a, double b) {
        return a <= b || aeq(a, b);
    }

    template <int dim=2>
    class Point<dim>;

    template<int dim>
    inline double dot(Point<dim> const& p1, Point<dim> const& p2) {
        double val = 0.0;
        for (size_t i = 0; i < dim + 1; i++) {
            val += p1[i] * p2[i];
        }
        return val;
    }

    inline double det2(double a1, double b1, double a2, double b2) {
        return a1 * b2 - b1 * a2;
    }


    inline Point<> intersection(Point<> const& p1, Point<> const& p2) {
        return Point<>(det2(p1[1], p1[2], p2[1], p2[2]),
                       det2(p1[0], p1[2], p2[0], p2[2]),
                       det2(p1[0], p1[1], p2[0], p2[1]));
    }

    template <int dim=2>
    class Point {
        // Points are of form (x, y, c) where c is the homgenous coordinate.
        std::array<double, dim + 1> coords;
    public:

        explicit template <typename ...Coords>
        Point(Coords... rest) : coords(rest...) {}

        Point() {
            coords.fill(0);
        }


        virtual std::string toString() const {
            std::ostringstream ss;
            ss << "pyscan::Point(";
            for (auto &el : coords) {
                ss << ", " << el;
            }
            ss << ")";
            return ss.str();
        }

        virtual bool operator==(Point<dim> const& pt) {

            for (int i = 0; i < dim + 1; i++) {
                if (coords[i] != pt.coords[i]) {
                    return false;
                }
            }
            return true;
        }

        virtual double& operater[](size_t i) {
            return &coords[i];
        }
        virtual double operator[](size_t i) const {
           return coords[i];
        }

        bool above_closed(Point<dim> const& p1) {
            return alte(dot(*this, p1), 0.0);
        }

        bool above(Point<dim> const& p1) {
            return alt(dot(*this, p1), 0.0);
        }

        bool parallel(Point<dim> const& l) {
            for (int i = 0; i < dim; i++) {
                if (!aeq(coords[i], l[i])) {
                    return false;
                }
            }
            return true;
        }
    };


    using point_list = std::vector<Point<>>;
    using weight_list = std::vector<double>;

}
#endif //PYSCAN_POINT_HPP
