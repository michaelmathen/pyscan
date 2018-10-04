//
// Created by mmath on 10/1/18.
//

#ifndef PYSCAN_RANGE_HPP
#define PYSCAN_RANGE_HPP

#include <unordered_set>
#include "Point.hpp"


namespace pyscan {


    template<int dim>
    class Range {
    public:
        virtual bool contains(Point <dim> const &pt) const = 0;

    };


    template<int dim>
    double range_weight(Range<dim> const& range, std::vector<WPoint<dim>> const& pts) {
        double weight = 0, total_weight = 0;
        for (auto const& pt : pts) {
            if (range.contains(pt)) {
                weight += pt.get_weight();
            }
            total_weight += pt.get_weight();
        }
        return weight / total_weight;
    }


    template<int dim>
    double range_weight(Range<dim> const& range, std::vector<LPoint<dim>> const& pts) {
        std::unordered_set<size_t> seen_labels;
        double weight = 0 , total_weight = 0;
        for (auto& pt : pts) {
            if ( (seen_labels.find(pt.get_label()) == seen_labels.end())) {
                weight += range.contains(pt) ? pt.get_weight() : 0.0;
                total_weight += pt.get_weight();
                seen_labels.emplace(pt.get_label());
            }
        }
        return weight;
    }

    template<typename R, typename Pt>
    double evaluate_range(R const& range,
                            std::vector<Pt> const& red_pts,
                            std::vector<Pt> const& blue_pts,
                            const discrepancy_func_t & scan) {
        return scan(range_weight(range, red_pts), range_weight(range, blue_pts));
    }


    template<typename R, typename NPt, typename SPt>
    auto scan_ranges2(std::vector<NPt> const &net,
                      std::vector<SPt> const &red_pts,
                      std::vector<SPt> const &blue_pts,
                      const discrepancy_func_t &scan) -> std::tuple<R, double> {

        double max_scan = -std::numeric_limits<double>::infinity();
        R max_region;
        for (auto p1 = net.begin(); p1 != net.end() - 1; p1++) {
            for (auto p2 = p1 + 1; p2 != net.end(); p2++) {
                R region(*p1, *p2);
                double curr_scan = evaluate_range(region, red_pts, blue_pts, scan);
                if (curr_scan > max_scan) {
                    max_region = region;
                    max_scan = curr_scan;
                }
            }
        }
        return std::make_tuple(max_region, max_scan);
    }

    /*
    * Will scan all combinatorial regions of type R that are defined by 3 points and return the region with the
    * highest discrepancy
    */
    template<typename R, typename Pt, typename SPt>
    auto scan_ranges3(std::vector<Pt> const &net,
                      std::vector<SPt> const &red_pts,
                      std::vector<SPt> const &blue_pts,
                      const discrepancy_func_t &scan) -> std::tuple<R, double> {
        double max_scan = -std::numeric_limits<double>::infinity();
        R max_region;
        for (auto p1 = net.begin(); p1 != net.end() - 2; p1++) {
            for (auto p2 = p1 + 1; p2 != net.end() - 1; p2++) {
                for (auto p3 = p2 + 1; p3 != net.end(); p3++) {
                    R region(*p1, *p2, *p3);
                    double curr_scan = evaluate_range(region, red_pts, blue_pts, scan);
                    if (curr_scan > max_scan) {
                        max_region = region;
                        max_scan = curr_scan;
                    }
                }
            }
        }
        return std::make_tuple(max_region, max_scan);
    }


    /*
     * Will scan all combinatorial regions of type R that are defined by 4 points and return the region with the
     * highest discrepancy
     */
    template<typename R, typename Pt, typename SPt>
    auto scan_ranges4(std::vector<Pt> const &net,
                      std::vector<SPt> const &red_pts,
                      std::vector<SPt> const &blue_pts,
                      const discrepancy_func_t &scan) -> std::tuple<R, double> {
        double max_scan = -std::numeric_limits<double>::infinity();
        R max_region;
        for (auto p1 = net.begin(); p1 != net.end() - 3; p1++) {
            for (auto p2 = p1 + 1; p2 != net.end() - 2; p2++) {
                for (auto p3 = p2 + 1; p3 != net.end() - 1; p3++) {
                    for (auto p4 = p3 + 1; p4 != net.end(); p4++) {
                        R region(*p1, *p2, *p3, *p4);
                        double curr_scan = evaluate_range(region, red_pts, blue_pts, scan);
                        if (curr_scan > max_scan) {
                            max_region = region;
                            max_scan = curr_scan;
                        }
                    }
                }
            }
        }
        return std::make_tuple(max_region, max_scan);
    }
}
#endif //PYSCAN_RANGE_HPP
