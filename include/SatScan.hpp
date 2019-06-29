/*
 * Created by Michael Matheny on 6/3/19.
 * at the University of Utah
 * email: michaelmathen@gmail.com
 * website: https://mmath.dev/
 */

#ifndef PYSCAN_SATSCAN_HPP
#define PYSCAN_SATSCAN_HPP
namespace pyscan {

    /*
     * These methods use disks that are defined by points in a grid or
     * as points in the measured and baseline sets.
     */
    std::tuple<Disk, double> satscan_grid(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double disk_r,
            discrepancy_func_t const& func);

    std::tuple<Disk, double> satscan_grid_labeled(
            const lpoint_list_t &measured,
            const lpoint_list_t &baseline,
            double grid_res,
            double disk_r,
            discrepancy_func_t const& func);

}
#endif //PYSCAN_SATSCAN_HPP
