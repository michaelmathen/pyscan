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


    template<typename T, typename Filter>
    double computeLabelTotalF(T begin, T end, std::unordered_map<size_t, size_t>& label_map, Filter filter) {
        double total = 0;
        for (; begin != end; ++begin) {
            if (filter(*begin)) {
                if (label_map.end() == label_map.find(begin->get_label())) {
                    total += begin->get_weight();
                    label_map[begin->get_label()] = 0;
                }
                label_map[begin->get_label()] = label_map[begin->get_label()] + 1;
            }
        }
        return total;
    }

    template< typename T>
    double computeLabelTotal(T begin, T end, std::unordered_map<size_t, size_t>& label_map) {
        return computeLabelTotalF(begin, end, label_map, [](pt2_t const& pt){ (void)pt; return true; });
    }

    template<typename T, typename Filter>
    double computeLabelTotalF(T begin, T end, Filter filter) {
        std::unordered_map<size_t, size_t> label_map;
        return computeLabelTotalF(begin, end, label_map, filter);
    }

    template< typename T>
    double computeLabelTotal(T begin, T end) {
        std::unordered_map<size_t, size_t> label_map;
        return computeLabelTotal(begin, end, label_map);
    }

    double updateCounts(std::unordered_map<size_t, size_t>& curr_counts,
                        crescent_t const& adding, crescent_t const& removing);

    template<typename T>
    double computeTotal(T begin, T end) {
        double sum = 0;
        std::for_each(begin, end, [&](WPoint<> const &pt) {
            sum += pt.get_weight();
        });
        return sum;
    }


    double evaluateRegion(
            wpoint_list_t const& m_pts,
            wpoint_list_t const& b_pts,
            Disk const& disk,
            const discrepancy_func_t&  scan);

    double evaluateRegion(
            lpoint_list_t const& m_pts,
            lpoint_list_t const& b_pts,
            Disk const& disk,
            const discrepancy_func_t& scan);

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


    std::tuple<Disk, double> diskScanSlow(
            point_list_t const& net,
            wpoint_list_t const& sampleM,
            wpoint_list_t const& sampleB,
            const discrepancy_func_t& scan);

    std::tuple<Disk, double> diskScanLabels(
            point_list_t net,
            lpoint_list_t sampleM,
            lpoint_list_t sampleB,
            const discrepancy_func_t& scan);



    std::tuple<Disk, double> diskScanSlowLabels(
            point_list_t net,
            lpoint_list_t sampleM,
            lpoint_list_t sampleB,
            const discrepancy_func_t& scan);


    std::tuple<Disk, double> diskScan(point_list_t net,
                  wpoint_list_t sampleM,
                  wpoint_list_t sampleB,
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


    std::tuple<Disk, double> disk_scan_simp(
            point_list_t const& net,
            wpoint_list_t const& sampleM,
            wpoint_list_t const& sampleB,
            const discrepancy_func_t& scan);


}
#endif //PYSCAN_DISKSCAN_HPP_HPP
