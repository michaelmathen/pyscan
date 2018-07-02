//
// Created by mmath on 5/28/17.
//

#ifndef PYSCAN_POINT_HPP
#define PYSCAN_POINT_HPP
#include <stdexcept>
#include <string>
#include <sstream>
#include <ostream>
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

    template <int dim> class Point;

    template<int dim> inline double dot(Point<dim> const&, Point<dim> const&);

    inline double det2(double a1, double b1, double a2, double b2);

    inline Point<2> intersection(Point<2> const& p1, Point<2> const& p2);

    template <int dim=2>
    class Point {
        // Points are of form (x, y, c) where c is the homgenous coordinate.
        std::array<double, dim + 1> coords;
    public:

        template <typename ...Coords>
        Point(Coords... rest) {
            std::initializer_list<double> il( {rest...} );
            std::copy (il.begin(), il.end(), coords.begin());
        }

        Point() {
            coords.fill(0);
        }

        virtual std::string str() const {
            std::stringstream ss;
            ss << *this;
            return ss.str();
        }

        double normalization_factor(Point<dim> const& pt) const {
            /*
             * Assumes that these are parallel. pt = (ax, ay, ac) so this returns
             * a
             */
            for (size_t i = 0; i < dim; i++) {
                if (!aeq(coords[i], 0.0)) {
                    return pt[i] / coords[i];
                }
            }
            return 1;
        }

        friend std::ostream& operator<<(std::ostream& os, Point const& pt)  { 
            os << "Point(";
            for (auto &el : pt.coords) {
                os << el << ", ";
            }
            os << ")";
            return os;
        }

        virtual bool operator==(Point<dim> const& pt) const {
            for (int i = 0; i < dim + 1; i++) {
                if (coords[i] != pt.coords[i]) {
                    return false;
                }
            }
            return true;
        }


        virtual bool approx_eq(Point<dim> const& pt) const {
            double scale1 = coords[dim];
            double scale2 = pt[dim];
            if (aeq(coords[dim], 0) && aeq(pt[dim], 0)) {
                scale1 = 1;
                scale2 = 1;
            } else if (aeq(coords[dim], 0) != aeq(pt[dim], 0)) {
                return false;
            }

            for (int d = 0; d < dim; d++) {
                if (!aeq(coords[d] / scale1, pt[d] / scale2)) {
                    return false;
                }
            }
            return true;
        }

        // virtual double& operator[](size_t i) {
        //     return &coords[i];
        // }
        virtual double operator[](size_t i) const {
           return coords[i];
        }

        virtual bool above_closed(Point<dim> const& p1) const {
            return alte(dot(*this, p1), 0.0);
        }

        virtual bool below_closed(Point<dim> const& p1) const {
            return alte(0.0, dot(*this, p1));
        }

        virtual bool above(Point<dim> const& p1) const {
            return alt(dot(*this, p1), 0.0);
        }

        virtual bool below(Point<dim> const& p1) const {
            return alt(0.0, dot(*this, p1));
        }


    };

    inline bool parallel(Point<2> const& l1, Point<2> const& l2) {
        return aeq(intersection(l1, l2)[2], 0);
    }

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


    inline Point<2> intersection(Point<2> const& p1, Point<2> const& p2) {
        return {det2(p1[1], p1[2], p2[1], p2[2]),
                       -det2(p1[0], p1[2], p2[0], p2[2]),
                       det2(p1[0], p1[1], p2[0], p2[1])};
    }


    inline double order_function(Point<> const& p1, Point<> const& p2) {
        double y = p1[2] * p2[0] - p2[2] * p1[0];
        double x = p1[2] * p2[1] - p2[2] * p1[1];
        if (y >= 0) {
            return atan2(y, x);
        } else {
            return M_PI + atan2(y, x);
        }
    }

    inline Point<> correct_orientation(Point<> const& pivot, Point<> const& p) {
        double y = pivot[2] * p[0] - p[2] * pivot[0];
        if (y >= 0) {
            return intersection(pivot, p);
        } else {
            return intersection(p, pivot);
        }
    }


    using Pt2 = Point<2>;

    using point_list = std::vector<Point<>>;
    using point_it = point_list::iterator;
    using weight_list = std::vector<double>;
    using weight_it = weight_list::iterator;
    using Pt3 = Point<3>;
    using point3_list = std::vector<Pt3>;
    using point3_it = point3_list::iterator;
    using label_list = std::vector<long>;


    inline void getLoc(Pt2 const& p, double& x, double& y) {
        x = p[0] / p[2];
        y = p[1] / p[2];
    }


    template<int dim>
    inline double getX(Point<dim> const& p1) {
        static_assert(dim > 0, "getX can only be used when dim > 0");
        return p1[0] / p1[dim];
    }

    template<int dim>
    inline double getY(Point<dim> const& p1) {
        static_assert(dim > 1, "getY can only be used when dim > 1");
        return p1[1] / p1[dim];
    }

    template<int dim>
    inline double getZ(Point<dim> const& p1) {
        static_assert(dim > 2, "getZ can only be used when dim > 2");
        return p1[2] / p1[dim];
    }

    template <int dim>
    inline bool sameLoc(Point<dim> const& p1, Point<dim> const& p2) {
        return p1.approx_eq(p2);
    }


    template<int dim=2>
    class WPoint : public Point<dim> {
    public:

        double weight;
        template <typename ...Coords>
        WPoint(double weight, Coords... rest) : Point<dim>(rest...), weight(weight) {
        }

        WPoint() : Point<dim>() {
            weight = 0;
        }
    };

    template<int dim=2>
    class LPoint : public WPoint<dim> {
    public:

        size_t label;
        template <typename ...Coords>
        LPoint(size_t label, double weight, Coords... rest) : WPoint<dim>(weight, rest...), label(label) {
        }

        LPoint() : WPoint<dim>() {
            label = 0;
        }
    };

    template<int dim>
    inline double getWeight(WPoint<dim> const& p) {
        return p.weight;
    }

    using wpoint_list = std::vector<WPoint<>>;
    using lpoint_list = std::vector<LPoint<>>;
}
#endif //PYSCAN_POINT_HPP
