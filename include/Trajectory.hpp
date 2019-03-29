//
// Created by mmath on 12/5/18.
//

#ifndef PYSCAN_TRAJECTORY_HPP
#define PYSCAN_TRAJECTORY_HPP

#include "Point.hpp"


namespace pyscan {

    class Trajectory {
        using pt_t = Point<2>;
        using p_list_t = std::vector<pt_t>;
        using p_it_t = typename p_list_t::iterator;
        using cp_it_t = typename p_list_t::const_iterator;

        p_list_t trajectory_pts;

    public:
        explicit Trajectory(point_list_t pts) : trajectory_pts(std::move(pts)) {}

        explicit Trajectory(std::vector<std::tuple<double, double>> const& pts) {
            trajectory_pts.reserve(pts.size());
            for (auto& p : pts) {
                trajectory_pts.emplace_back(p);
            }
        }


        virtual double get_length() const {
            if (trajectory_pts.empty()) return 0.0;

            double total_distance = 0;
            auto st_pt = trajectory_pts[0];
            for (size_t i = 1; i < trajectory_pts.size(); i++) {
                total_distance += st_pt.dist(trajectory_pts[i]);
                st_pt = trajectory_pts[i];
            }
            return total_distance;
        }

        virtual ~Trajectory() = default;

        virtual double get_weight() const {
            return 1;
        }

        p_list_t get_pts() const {
            return trajectory_pts;
        }

        size_t size() const {
            return trajectory_pts.size();
        }

        const pt_t &operator[](size_t i) const {
            return trajectory_pts[i];
        }


        virtual double get_partial_weight() const {
            return get_length();
        }

        cp_it_t begin() const {
            return trajectory_pts.begin();
        }

        cp_it_t end() const {
            return trajectory_pts.end();
        }

        p_it_t begin() {
            return trajectory_pts.begin();
        }

        p_it_t end() {
            return trajectory_pts.end();
        }

        bool empty() const {
            return trajectory_pts.empty();
        }

        virtual bool operator==(const Trajectory& other) {
            if (other.size() == this->size() || other.get_weight() != this->get_weight()) {
                for (size_t i = 0; i < other.size(); i++) {
                    if (!(other.trajectory_pts[i] == trajectory_pts[i])) {
                        return false;
                    }
                }
                return true;
            } else {
                return false;
            }
        }

        inline double point_dist(const Point<2> &p1) const {
            /*
             * Computes the distance to the trajectory.
             */
            if (trajectory_pts.empty()) {
                return std::numeric_limits<double>::infinity();
            } else if (trajectory_pts.size() == 1) {
                return p1.dist(trajectory_pts[0]);
            } else {
                double min_dist = std::numeric_limits<double>::infinity();
                auto last_pt = trajectory_pts.begin();
                for (auto curr_pt = last_pt + 1; curr_pt != trajectory_pts.end(); ++curr_pt) {
                    double seg_dist = p1.square_dist(*last_pt, *curr_pt);
                    if (seg_dist < min_dist) {
                        min_dist = seg_dist;
                    }
                    last_pt = curr_pt;
                }
                return sqrt(min_dist);
            }
        }
    };

    class WTrajectory : public Trajectory {
        double weight;
    public:
        WTrajectory(double w, point_list_t pts) : Trajectory(std::move(pts)), weight(w) {}
        WTrajectory(double w, std::vector<std::tuple<double, double>> const& pts): Trajectory(pts), weight(w) {}

        double get_weight() const override {
            return weight;
        }

        double get_partial_weight() const override {
            return weight * this->get_length();
        }

    };

    using trajectory_t = Trajectory;
    using wtrajectory_t = WTrajectory;
    using trajectory_set_t = std::vector<trajectory_t>;
    using wtrajectory_set_t = std::vector<wtrajectory_t>;
}
#endif //PYSCAN_TRAJECTORY_HPP
