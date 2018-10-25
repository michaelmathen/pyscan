#ifndef __RANGE_H__
#define __RANGE_H__

#include "Common.hpp"
#include <unordered_set>

namespace pyscan {

template <int dim>
class Range {
public:
    virtual bool contains(const Point<dim>& pt) const = 0;
    virtual bool intersects_segment(const Point<dim> &p1, const Point<dim> &p2) const = 0;
};

template <int dim>
inline double computeTotal(const std::vector<WPoint < dim>>& pts) {
    double res = 0.0;
    for (auto& x: pts) res += x.get_weight();
    return res;
}

template <int dim>
inline double computeTotal(const std::vector<LPoint < dim>>& pts) {
    double res = 0.0;
    std::unordered_set<size_t> seen;
    for (auto& x: pts) {
        if (seen.find(x.get_label()) == seen.end()) {
            res += x.get_weight();
            seen.emplace(x.get_label());
        }
    }
    return res;
}

template <int dim>
double range_weight(const Range<dim>& range, const std::vector<WPoint<dim>>& pts) {
    double weight = 0.0;
    for (auto& pt: pts) {
        if (range.contains(pt)) {
            weight += pt.get_weight();
        }
    }
    return weight;
}

template <int dim>
double range_weight(const Range<dim>& range, const std::vector<LPoint<dim>>& pts) {
    std::unordered_set<size_t> seen;
    double weight = 0.0;
    for (auto& pt: pts) {
        if (range.contains(pt) && seen.find(pt.get_label()) == seen.end()) {
            weight += pt.get_weight();
            seen.emplace(pt.get_label());
        }
    }
    return weight;
}

template <int dim, typename Pt>
double evaluate_range(
        const Range<dim>& range,
        const std::vector<Pt>& red,
        const std::vector<Pt>& blue,
        const discrepancy_func_t& f) {

    return f(range_weight(range, red), computeTotal(red),
            range_weight(range, blue), computeTotal(blue));
}

template <typename R, int dim = 2>
std::tuple<R, double> max_range2(
        const point_list_t &point_net,
        const std::vector<WPoint < dim>>& red,
        const std::vector<WPoint<dim>>& blue,
        const discrepancy_func_t& f) {

    double max_stat = 0.0;
    R cur_max;
    for (size_t i = 0; i < point_net.size() - 2; ++i) {
        for (size_t j = i + 1; j < point_net.size() - 1; ++j) {
            R now(point_net[i], point_net[j]);
            double cur_stat = evaluate_range(now, red, blue, f);
            if (cur_stat > max_stat) {
                cur_max = now;
                max_stat = cur_stat;
            }
        }
    }

    return std::make_tuple(cur_max, max_stat);
}

template <typename R, int dim = 2>
std::tuple<R, double> max_range2_labeled(
        const point_list_t &point_net,
        const std::vector<LPoint < dim>>& red,
        const std::vector<LPoint<dim>>& blue,

        const discrepancy_func_t& f) {

    double max_stat = 0.0;
    R cur_max;
    for (size_t i = 0; i < point_net.size() - 2; ++i) {
        for (size_t j = i + 1; j < point_net.size() - 1; ++j) {
            R now(point_net[i], point_net[j]);
            double cur_stat = evaluate_range(now, red, blue, f);
            if (cur_stat > max_stat) {
                cur_max = now;
                max_stat = cur_stat;
            }
        }
    }
    return std::make_tuple(cur_max, max_stat);
}

template <typename R, int dim = 2>
std::tuple<R, double> max_range3(
        const std::vector<Point<dim>> &point_net,
        const std::vector<WPoint < dim>>& red,
        const std::vector<WPoint<dim>>& blue,
        const discrepancy_func_t& f) {

    double max_stat = 0.0;
    R cur_max;
    for (size_t i = 0; i < point_net.size() - 2; ++i) {
        for (size_t j = i + 1; j < point_net.size() - 1; ++j) {
            for (size_t k = j + 1; k < point_net.size(); ++k) {
                R now(point_net[i], point_net[j], point_net[k]);
                double cur_stat = evaluate_range(now, red, blue, f);
                if (cur_stat > max_stat) {
                    cur_max = now;
                    max_stat = cur_stat;
                }
            }
        }
    }

    return std::make_tuple(cur_max, max_stat);
}

template <typename R, int dim = 2>
std::tuple<R, double> max_range3_labeled(
        const std::vector<Point<dim>> &point_net,
        const std::vector<LPoint < dim>>& red,
        const std::vector<LPoint<dim>>& blue,
        const discrepancy_func_t& f) {

    double max_stat = 0.0;
    R cur_max;
    for (size_t i = 0; i < point_net.size() - 2; ++i) {
        for (size_t j = i + 1; j < point_net.size() - 1; ++j) {
            for (size_t k = j + 1; k < point_net.size(); ++k) {
                R now(point_net[i], point_net[j], point_net[k]);
                double cur_stat = evaluate_range(now, red, blue, f);
                if (cur_stat > max_stat) {
                    cur_max = now;
                    max_stat = cur_stat;
                }
            }
        }
    }

    return std::make_tuple(cur_max, max_stat);
}

}

#endif
