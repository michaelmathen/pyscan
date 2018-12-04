#ifndef __DISK_SCAN_H__
#define __DISK_SCAN_H__

#include "Common.hpp"
#include "Point.hpp"
#include "Disk.hpp"

namespace pyscan {

    std::tuple<Disk, double> max_disk(
            const point_list_t &point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            const discrepancy_func_t &f);

    std::tuple<Disk, double> max_disk_labeled(
            const point_list_t &point_net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            const discrepancy_func_t &f);

//    std::tuple<Disk, double> max_disk_lift(
//            const point_list_t &point_net,
//            const wpoint_list_t &red,
//            const wpoint_list_t &blue,
//            const discrepancy_func_t &f);
//
//    std::tuple<Disk, double> max_disk_lift_labeled(
//            const point_list_t &point_net,
//            const lpoint_list_t &red,
//            const lpoint_list_t &blue,
//            const discrepancy_func_t &f);


    std::tuple<Disk, double> max_disk_scale(
            const point_list_t &point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            double min_res,
            const discrepancy_func_t &f);

    std::tuple<Disk, double> max_disk_scale_labeled(
            const lpoint_list_t &point_net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            bool compress,
            double min_res,
            const discrepancy_func_t &f);

    std::tuple<Disk, double> max_disk_scale_labeled_alt(
            const point_list_t &point_net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            double min_res,
            const discrepancy_func_t &f);

    std::tuple<Disk, double> max_disk_scale_slow(
            const point_list_t &point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            double min_res,
            double max_res,
            const discrepancy_func_t &f);

    std::tuple<Disk, double> max_disk_scale_slow_labeled(
            const point_list_t &point_net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            double min_res,
            double max_res,
            const discrepancy_func_t &f);


#ifdef DEBUG

    //Uses an internal
    std::tuple<Disk, double> max_disk_scale_alt(
            const lpoint_list_t &point_net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            double min_res,
            double max_res,
            const discrepancy_func_t &f);
#endif

    //    std::tuple<Disk, double> max_rdisk_lift_labeled(
//            const lpoint_list_t &net,
//            const lpoint_list_t &red,
//            const lpoint_list_t &blue,
//            double alpha,
//            double min_radius,
//            double max_radius,
//            const discrepancy_func_t &f);
//}

}
#endif
