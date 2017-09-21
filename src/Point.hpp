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

    template <typename Weight=int>
    class Point {
        double x;
        double y;
        Weight red;
        Weight blue;
    public:

        Point(double x, double y, Weight r, Weight b) : x(x), y(y), red(r), blue(b) {}

        Weight getWeight() const {
            return red + blue;
        }
        Weight getRedWeight() const {
            return red;
        }
        Weight getBlueWeight() const {
            return blue;
        }
        double getX() const {
            return x;
        }
        double getY() const {
            return y;
        }

        std::string toString() const {
            std::ostringstream ss;
            ss << "pyscan::Point(" <<  x << ", " << y;
            ss << ", " << red << ", " << blue;
            ss << ")";
            return ss.str();
        }

    };


    using point_it = std::vector<Point<int> >::iterator;
    using point_d_it = std::vector<Point<double>>::iterator;

}
#endif //PYSCAN_POINT_HPP
