//
// Created by mmath on 9/4/18.
//


#ifndef PYSCAN_TRAJECTORYSCAN_HPP
#define PYSCAN_TRAJECTORYSCAN_HPP

#include "RectangleScan.hpp"
#include "Disk.hpp"
#include "HalfSpaceScan.hpp"
#include "Point.hpp"

namespace pyscan {



    class Trajectory {
        point_list_t trajectory_pts;
    public:
        explicit Trajectory(point_list_t pts) : trajectory_pts(std::move(pts)) {}

        double get_length() const {
            if (trajectory_pts.empty()) return 0.0;

            double total_distance = 0;
            auto st_pt = trajectory_pts[0];
            for (size_t i = 1; i < trajectory_pts.size(); i++) {
                total_distance += st_pt.dist(trajectory_pts[i]);
                st_pt = trajectory_pts[i];
            }
            return total_distance;
        }

        virtual double get_weight() const {
            return 1;
        }

        point_list_t get_pts () const {
            return trajectory_pts;
        }
        
        virtual double get_partial_weight() const {
            return get_length();
        }

        cpoint_it_t begin() const {
            return trajectory_pts.begin();
        }

        cpoint_it_t end() const {
            return trajectory_pts.end();
        }

        point_it_t begin() {
            return trajectory_pts.begin();
        }

        point_it_t end() {
            return trajectory_pts.end();
        }

        bool empty() const {
            return trajectory_pts.empty();
        }

        inline double point_dist(const pt2_t &p1) const {
            /*
             * Computes the distance to the trajectory.
             */
            if (trajectory_pts.empty()) {
                return std::numeric_limits<double>::infinity();
            } else if (trajectory_pts.size() == 1){
                return p1.dist(trajectory_pts[0]);
            } else {
                double min_dist = std::numeric_limits<double>::infinity();
                auto last_pt = trajectory_pts.begin();
                for (auto curr_pt = last_pt + 1; curr_pt != trajectory_pts.end(); ++curr_pt) {
                    double seg_dist = p1.square_dist(*last_pt, *curr_pt);
                    if (seg_dist < min_dist) {
                        min_dist = seg_dist;
                    }
                }
                return min_dist;
            }
        }


        inline bool intersects_region(const Range<2>& range) {
            if (trajectory_pts.empty()) {
                return false;
            } else if (trajectory_pts.size() == 1){
                return range.contains(trajectory_pts[0]);
            } else {
                auto last_pt = trajectory_pts.begin();
                for (auto curr_pt = last_pt + 1; curr_pt != trajectory_pts.end(); ++curr_pt) {
                    if (range.intersects_segment(*last_pt, *curr_pt)) return true;
                }
                return false;
            }
        }

        inline bool intersects_disk(const Disk& range) {
            return intersects_region(range);
        }

        inline bool intersects_halfplane(const HalfSpace<2>& range) {
            return intersects_region(range);
        }
    };

    class WTrajectory : public Trajectory {
        double weight;
    public:
        WTrajectory(double w, point_list_t pts) : Trajectory(std::move(pts)), weight(w) {}

        double get_weight() const override {
            return weight;
        }

        double get_partial_weight() const override {
            return weight * get_length();
        }

    };

    using trajectory_t = Trajectory;
    using wtrajectory_t = WTrajectory;
    using trajectory_set_t = std::vector<trajectory_t>;
    using wtrajectory_set_t = std::vector<wtrajectory_t>;


    //////////////////////////////////////
    //Full Scanning code//////////////////
    //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Disk scanning Trajectory code//////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::tuple<Disk, double> max_disk_traj_grid(trajectory_set_t const& net,
                                                wtrajectory_set_t const& sampleM,
                                                wtrajectory_set_t const& sampleB,
                                                double alpha,
                                                double max_r,
                                                const discrepancy_func_t &scan);

}
    


#endif //PYSCAN_TRAJECTORYSCAN_HPP
