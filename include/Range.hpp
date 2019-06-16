#ifndef __RANGE_H__
#define __RANGE_H__

#include "Common.hpp"
#include "Trajectory.hpp"

#include <unordered_set>

namespace pyscan {

template <int dim>
class Range {
public:
    virtual bool contains(const Point<dim>& pt) const = 0;
    virtual bool intersects_segment(const Point<dim> &p1, const Point<dim> &p2) const = 0;

    template <typename Dummy = bool>
    typename std::enable_if<dim == 2, Dummy>::type intersects_trajectory(Trajectory const& trajectory) const {
        if (trajectory.empty()) {
            return false;
        } else if (trajectory.size() == 1) {
            return contains(trajectory[0]);
        } else {
            auto last_pt = trajectory.begin();
            for (auto curr_pt = last_pt + 1; curr_pt != trajectory.end(); ++curr_pt) {
                if (intersects_segment(*last_pt, *curr_pt)) {
                    return true;
                }
                last_pt = curr_pt;
            }
            return false;
        }
    }
    virtual ~Range() = default;
};

template<typename Traj=WTrajectory>
double computeTotal(const std::vector<Traj> &traj_set) {
    double weight = 0;
    for (auto& traj : traj_set) {
        weight += traj.get_weight();
    }
    return weight;
}

template <typename Traj=WTrajectory>
typename std::enable_if<std::is_base_of<Trajectory, Traj>::value, double>::type
range_weight(const Range<2>& range,
                    const std::vector<Trajectory>& traj_set) {
    double weight = 0;
    for (auto& traj : traj_set) {
        if (range.intersects_trajectory(traj)) {
            weight += traj.get_weight();
        }
    }
    return weight;
}

template <int dim>
inline double computeTotal(const std::vector<WPoint < dim>>& pts) {
    double res = 0.0;
    for (auto& x: pts) res += x.get_weight();
    return res;
}

template <typename It>
inline double computeTotal(It b, It e) {
    double res = 0.0;
    for (;b!=e; b++) res += b->get_weight();
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

template <int dim, typename Obj>
double evaluate_range(
        const Range<dim>& range,
        const std::vector<Obj>& red,
        const std::vector<Obj>& blue,
        const discrepancy_func_t& f) {

    return f(range_weight(range, red), computeTotal(red),
            range_weight(range, blue), computeTotal(blue));
}

template <typename R, template <int> typename P=WPoint, int dim=2>
std::tuple<R, double> max_range2(
        const std::vector<Point<dim>> &point_net,
        const std::vector<P<dim>>& red,
        const std::vector<P<dim>>& blue,
        const discrepancy_func_t& f) {

    double max_stat = 0.0;
    R cur_max;
    for (size_t i = 0; i < point_net.size() - 1; ++i) {
        for (size_t j = i + 1; j < point_net.size(); ++j) {
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

template <typename R, template <int> typename P=WPoint, int dim=2>
std::tuple<R, double> max_range3(
        const std::vector<Point<dim>>& point_net,
        const std::vector<P<dim>>& red,
        const std::vector<P<dim>>& blue,
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
