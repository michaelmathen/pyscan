//
// Created by mmath on 5/28/17.
//

#ifndef PYSCAN_POINT_HPP
#define PYSCAN_POINT_HPP
#include <stdexcept>
#include <initializer_list>
#include <string>
#include <sstream>
#include <vector>

#include "Vecky.hpp"

namespace pyscan {

    template <int dim=2>
    class Point {
        double red;
        double blue;

        VecN<double, dim> coords;
    public:

        template <typename ...Coords>
        Point(double r, double b, Coords... rest) : coords(rest...) {
            red = r;
            blue = b;
            static_assert(dim == sizeof...(Coords), "coords has to be the same as the number of dimensions");
        }

        Point() : red(0), blue(0){
        }

      	virtual void setRedWeight(double w_r) {
      	    red = w_r;
      	}

      	virtual void setBlueWeight(double w_b) {
      	    blue = w_b;
      	}

        virtual double getWeight() const {
            return red + blue;
        }
        virtual double getRedWeight() const {
            return red;
        }
        virtual double getBlueWeight() const {
            return blue;
        }

        virtual std::string toString() const {
            std::ostringstream ss;
            ss << "pyscan::Point(";
            ss << red << ", " << blue;
            for (int i = 0; i < dim; i++) {
                ss << ", " << coords[i];
            }
            ss << ")";
            return ss.str();
        }

        template <int ix,  int d>
        friend double get(Point<d> const& pt);

        virtual bool operator==(Point<dim> const& pt) {

            bool val = pt.getRedWeight() == getRedWeight() &&
                       pt.getBlueWeight() == getBlueWeight();
            if (!val)
                return false;
            for (int i = 0; i < dim; i++) {
                if (coords[i] != pt.coords[i]) {
                    return false;
                }
            }
            return true;
        }

        virtual double getX() const {
            return coords[0];
        }

        virtual double getY() const {
            return coords[1];
        }
        virtual double get(int i) const {
          return coords[i];
        }
    };


    template<int dim>
    inline double dot(Point<dim> const& p1, Point<dim> const& p2) {
        return dot(p1.coords, p2.coords);
    }

    template <int dim=2>
    class LPoint : public Point<dim> {
        size_t label;
    public:
        template <typename ...Coords>
        LPoint(size_t l, double r, double b, Coords... rest) : label(l), Point<dim>(r, b, rest...) {}

        virtual size_t getLabel() const {
            return label;
        }
        virtual bool operator==(LPoint<dim> const& lpt) {
            return Point<dim>::operator==(lpt) && lpt.label == this->label;
        }
    };

    template <int ix, int dim>
    double get(Point<dim> const& pt) {
        static_assert(ix < dim, "Requested a coordinate that is greater than the dimension");
        static_assert(0 <= ix, "Requested a coordinate that is less than zero");
        return pt.get(ix);
    }

    using point_it = std::vector<Point<>>::iterator;

    double getMeasured(Point<> const& pt);
    double getBaseline(Point<> const& pt);

    void getLoc(Point<> const& pt, double& x1, double& x2);

    double getX(Point<> const& pt);
    double getY(Point<> const& pt);
    bool colinear(Point<> const&, Point<> const&, Point<> const&);
    bool onLineSegment(Point<> const& pt1,
                       Point<> const& pt2,
                       Point<> const& pt3);

    bool sameLoc(Point<> const& p1, Point<> const& p2);


    inline Point<> removeLabel(LPoint<> const& pt) {
      return Point<>(pt.getRedWeight(), pt.getBlueWeight(), pt.get(0), pt.get(1));
    }

}
#endif //PYSCAN_POINT_HPP
