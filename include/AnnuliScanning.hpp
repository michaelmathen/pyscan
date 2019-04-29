/*
 * Created by Michael Matheny on 4/25/19.
 * at the University of Utah
 * email: michaelmathen@gmail.com
 * website: https://mmath.dev/
 */

#ifndef PYSCAN_ANNULISCANNING_HPP
#define PYSCAN_ANNULISCANNING_HPP

#include "Disk.hpp"
#include "Point.hpp"

namespace pyscan {
    std::tuple<Disk, double> max_annuli(const point_list_t &pts,
                                        const wpoint_list_t &mpts,
                                        const wpoint_list_t &bpts,
                                        const std::vector<double> &radii,
                                        const discrepancy_func_t &func);
}
#endif //PYSCAN_ANNULISCANNING_HPP
