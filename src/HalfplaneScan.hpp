//
// Created by mmath on 9/19/17.
//

#ifndef PYSCAN_HALFPLANE_HPP
#define PYSCAN_HALFPLANE_HPP

#include <functional>
#include <tuple>
#include "Point.hpp"

namespace pyscan {
    Halfspace<3> liftHalfspace(Halfspace<2> const& h2, Point<3> const& p3);

    Point<2> dropPoint(Point<3> const& fixed_point, Point<3> const& p3);

    double angle(Point<> const & p1, Point<> const& p2);

    std::tuple<Point<>, double> max_halfplane(point_list& red, weight_list& red_weight,
                                              point_list& blue, weight_list& blue_weight,
                                              std::function<double(double, double)> f);
}
#endif //PYSCAN_HALFPLANE_HPP
