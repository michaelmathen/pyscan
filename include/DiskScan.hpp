//
// Created by mmath on 9/25/17.
//

#ifndef PYSCAN_DISKSCAN_HPP_HPP
#define PYSCAN_DISKSCAN_HPP_HPP

#include <vector>
#include <unordered_map>
#include <algorithm>

#include "DiskScan2.hpp"
#include "Point.hpp"

namespace pyscan {


    struct LabeledVal {
        size_t label;
        double value;
        LabeledVal(size_t l, double v) : label(l), value(v) {}
    };
    using crescent_t = std::vector<LabeledVal>;



    double updateCounts(std::unordered_map<size_t, size_t>& curr_counts,
                        crescent_t const& adding, crescent_t const& removing);

    void solveCircle3(
            pt2_t const& pt1,
            pt2_t const& pt2,
            pt2_t const& pt3,
            double &a, double &b);

    void solveCircle3(
            pt2_t const& pt1,
            pt2_t const& pt2,
            pt2_t const& pt3,
            double &a, double &b, double &r);



    std::tuple<Disk, double> disk_scan(point_list_t const& net,
                  wpoint_list_t const& sampleM,
                  wpoint_list_t const& sampleB,
                  const discrepancy_func_t& scan);

    std::tuple<Disk, double> disk_scan_labels(point_list_t const& net,
                                       lpoint_list_t const& sampleM,
                                       lpoint_list_t const& sampleB,
                                       const discrepancy_func_t& scan);

    std::tuple<Disk, double> disk_scan_restricted(
            pt2_t p1,
            pt2_t p2,
            point_list_t net,
            wpoint_list_t const& sampleM,
            wpoint_list_t const& sampleB,
            double min_dist,
            double max_dist,
            double m_Total,
            double b_Total,
            const discrepancy_func_t& scan);

    std::tuple<Disk, double> disk_scan_restricted(
            pt2_t p1,
            pt2_t p2,
            point_list_t net,
            lpoint_list_t const& sampleM,
            lpoint_list_t const& sampleB,
            double min_dist,
            double max_dist,
            double m_Total,
            double b_Total,
            const discrepancy_func_t& scan);


    template<typename T>
    std::tuple<Disk, double> disk_scan_scale(
            point_list_t const& net,
            std::vector<T> const& sampleM,
            std::vector<T> const& sampleB,
            uint32_t grid_r,
            const discrepancy_func_t& scan);

    template<typename T>
    std::tuple<Disk, double> cached_disk_scan(
            point_list_t const& net,
            T const& sampleM,
            T const& sampleB,
            const discrepancy_func_t& scan);


    std::tuple<Disk, double> disk_scan_simple(
            point_list_t const& net,
            wpoint_list_t const& sampleM,
            wpoint_list_t const& sampleB,
            const discrepancy_func_t& scan);

    std::tuple<Disk, double> disk_scan_simple_labels(
            point_list_t const& net,
            lpoint_list_t const& sampleM,
            lpoint_list_t const& sampleB,
            discrepancy_func_t const& scan);
}
#endif //PYSCAN_DISKSCAN_HPP_HPP
