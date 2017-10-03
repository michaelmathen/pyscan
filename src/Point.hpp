//
// Created by mmath on 5/28/17.
//

#ifndef PYSCAN_POINT_HPP
#define PYSCAN_POINT_HPP
#include <stdexcept>
#include <initializer_list>
#include <string>
#include <sstream>

#include "Vecky.hpp"

namespace pyscan {

    template <typename Weight=int, int dim=2>
    class Point {
        Weight red;
        Weight blue;

        VecN<Weight, dim> coords;
    public:

        template <typename ...Coords>
        Point(Weight r, Weight b, Coords... rest) : coords(rest...) {
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

        virtual Weight getWeight() const {
            return red + blue;
        }
        virtual Weight getRedWeight() const {
            return red;
        }
        virtual Weight getBlueWeight() const {
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

        template <int ix, typename W, int d>
        friend double get(Point<W, d> const& pt);

        virtual bool operator==(Point<Weight, dim> const& pt) {

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
    };


    template<typename W, int dim>
    inline double dot(Point<W, dim> const& p1, Point<W, dim> const& p2) {
        return dot(p1.coords, p2.coords);
    }

    template <typename Weight=int, int dim=2>
    class LPoint : public Point<Weight, dim> {
        size_t label;
    public:
        template <typename ...Coords>
        LPoint(size_t label, Weight r, Weight b, Coords... rest) : Point<Weight, dim>(r, b, rest...) {}

        size_t getLabel() const {
            return label;
        }
        virtual bool operator==(LPoint<Weight, dim> const& lpt) {
            auto& pt1 = (Point<Weight, dim>&)lpt;
            auto& pt2 = (Point<Weight, dim>&)*this;
            return pt2 == pt2 && lpt.label == this->label;
        }
    };

    template <int ix, typename W, int dim>
    double get(Point<W, dim> const& pt) {
        static_assert(ix < dim, "Requested a coordinate that is greater than the dimension");
        return pt.coords[ix];
    }

    using point_it = std::vector<Point<int, 2>>::iterator;
    using point_d_it = std::vector<Point<double, 2>>::iterator;

}
#endif //PYSCAN_POINT_HPP
