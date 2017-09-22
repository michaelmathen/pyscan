//
// Created by mmath on 5/28/17.
//

#ifndef PYSCAN_POINT_HPP
#define PYSCAN_POINT_HPP
#include <stdexcept>
#include <initializer_list>
#include <string>
#include <sstream>

namespace pyscan {

    template <typename Weight=int, int dim=2>
    class Point {
        Weight red;
        Weight blue;

        double coords[dim];

        template <typename F>
        Point(int ix, F el) {
            coords[dim - ix] = el;
        }

        template <typename F, typename ...Coords>
        Point(int ix, F el, Coords... rest) : Point(ix - 1, rest...){
            coords[dim - ix] = el;
            static_assert(std::is_same<double, F>::value, "One of the coords is not of type double");
        }

    public:

        template <typename ...Coords>
        Point(Weight r, Weight b, Coords... rest) : Point(dim, rest...){
            red = r;
            blue = b;
            static_assert(dim == sizeof...(Coords), "coords has to be the same as the number of dimensions");
        }


        Weight getWeight() const {
            return red + blue;
        }
        Weight getRedWeight() const {
            return red;
        }
        Weight getBlueWeight() const {
            return blue;
        }

        template<int ix>
        double get() const {
            static_assert(ix < dim, "Requested a coordinate that is greater than the dimension");
            return coords[ix];
        }

        std::string toString() const {
            std::ostringstream ss;
            ss << "pyscan::Point(";
            ss << red << ", " << blue;
            for (int i = 0; i < dim; i++) {
                ss << ", " << coords[i];
            }
            ss << ")";
            return ss.str();
        }

    };


    using point_it = std::vector<Point<int, 2>>::iterator;
    using point_d_it = std::vector<Point<double, 2>>::iterator;

}
#endif //PYSCAN_POINT_HPP
