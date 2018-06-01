//
// Created by mmath on 5/22/18.
//

#ifndef PYSCAN_DISKSCAN2_HPP
#define PYSCAN_DISKSCAN2_HPP

#include "Point.hpp"

namespace pyscan {

    std::tuple<Pt3, double> max_disk(
            point3_list& point_net,
            point3_list& red,
            weight_list& red_w,
            point3_list& blue,
            weight_list& blue_w,
            std::function<double(double, double)> const& f);
}
#endif //PYSCAN_DISKSCAN2_HPP
