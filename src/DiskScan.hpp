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
                if (label_map.end() == label_map.find(begin->label)) {
                    total += getWeight(*begin);
                    label_map[begin->label] = 0;
                }
                label_map[begin->label] = label_map[begin->label] + 1;
            }
        }
        return total;
    }

    template< typename T>
    double computeLabelTotal(T begin, T end, std::unordered_map<size_t, size_t>& label_map) {
        return computeLabelTotalF(begin, end, label_map, [](Point<> const& pt){ return true; });
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
            sum += getWeight(pt);
        });
        return sum;
    }


    auto evaluateRegion(pyscan::wpoint_list const& m_pts, pyscan::wpoint_list const& b_pts,
                        pyscan::Disk const& disk, std::function<double(double, double)> const& scan)-> double;

    auto evaluateRegion(pyscan::lpoint_list const& m_pts, pyscan::lpoint_list const& b_pts,
                        pyscan::Disk const& disk, std::function<double(double, double)> const& scan)-> double;

    void solveCircle3(Point<> const& pt1, Point<> const& pt2, Point<> const& pt3,
                      double &a, double &b);

    void solveCircle3(Point<> const& pt1, Point<> const& pt2, Point<> const& pt3,
                      double &a, double &b, double &r);


    std::tuple<Disk, double> diskScanSlow(point_list const& net, wpoint_list const& sampleM, wpoint_list const& sampleB,
                                          std::function<double(double, double)> const& scan);

    std::tuple<Disk, double> diskScanLabels(point_list& net,
                        lpoint_list& sampleM,
                        lpoint_list& sampleB,
                        std::function<double(double, double)> const& scan);



    std::tuple<Disk, double> diskScanSlowLabels(point_list& net, lpoint_list& sampleM, lpoint_list& sampleB,
                            std::function<double(double, double)> const& scan);


    std::tuple<Disk, double> diskScan(point_list& net,
                  wpoint_list& sampleM,
                  wpoint_list& sampleB,
                  std::function<double(double, double)> const& scan);


    auto disk_scan_restricted(Point<> p1,
                              Point<>  p2,
                              std::vector<Point<2>>& net,
                              std::vector<WPoint<2>>& sampleM,
                              std::vector<WPoint<2>>& sampleB,
                              double min_dist,
                              double max_dist,
                              double m_Total,
                              double b_Total,
                              std::function<double(double, double)> const& scan)
    -> std::tuple<Disk, double>;

    template<typename T>
    std::tuple<Disk, double> disk_scan_scale(point_list& net,
                                        std::vector<T>& sampleM,
                                        std::vector<T>& sampleB,
                                        uint32_t grid_r,
                                        std::function<double(double, double)> const& scan);

    template<typename T>
    std::tuple<Disk, double> cached_disk_scan(point_list& net,
                                              T& sampleM,
                                              T& sampleB,
                                              std::function<double(double, double)> const& scan);



    std::tuple<Disk, double> disk_scan_simp(point_list& net,
                                            std::vector<WPoint<>>& sampleM,
                                            std::vector<WPoint<>>& sampleB,
                                            std::function<double(double, double)> const& scan);


}
#endif //PYSCAN_DISKSCAN_HPP_HPP
