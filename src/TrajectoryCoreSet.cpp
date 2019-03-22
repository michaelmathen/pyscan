//
// Created by mmath on 9/28/18.
//

#include <random>
#include <unordered_map>
#include <unordered_set>
#include <cassert>

#include "appext.h"
#include "Point.hpp"
#include "Disk.hpp"
#include "FunctionApprox.hpp"

#include "TrajectoryCoreSet.hpp"

namespace pyscan {




    /*
  * Compute the lower leftmost corner of a box containing these points.
  */
    std::tuple<double, double, double, double> bounding_box(point_list_t::const_iterator traj_b,
                                                            point_list_t::const_iterator traj_e) {

        double lx = (*traj_b)(0), ux = (*traj_b)(0), ly = (*traj_b)(1), uy = (*traj_b)(1);
        for ( ; traj_b != traj_e; traj_b++) {
            lx = std::min(lx, (*traj_b)(0));
            ux = std::max(ux, (*traj_b)(0));
            ly = std::min(ly, (*traj_b)(1));
            uy = std::max(uy, (*traj_b)(1));
        }
        return std::make_tuple(lx, ly, ux, uy);
    }


    double x_to_y(double x_1, double y_1, double x_2, double y_2, double x) {

        return (y_1 - y_2) / (x_1 - x_2) * (x - x_1) + y_1;
    }

    double y_to_x(double x_1, double y_1, double x_2, double y_2, double y) {
        return (x_1 - x_2) / (y_1 - y_2) * (y - y_1) + x_1;

    }

    long index(double x, double y, double lx, double ly, double chord_l, long g_size) {
        long i = static_cast<long>((x - lx) / chord_l);
        long j = static_cast<long>((y - ly) / chord_l);
        return i + g_size * j;
    }

    bool inside_box(double x, double y, pt2_t const& lpt, pt2_t const& rpt) {
        double ux = std::max(lpt(0), rpt(0));
        double lx = std::min(lpt(0), rpt(0));
        double uy = std::max(lpt(1), rpt(1));
        double ly = std::min(lpt(1), rpt(1));
        return util::alte(x, ux) && util::alte(lx, x) &&  util::alte(y, uy) && util::alte(ly, y);
    }
    /*
     * Takes a trajectory and grids it so that each grid contains points that cross the boundaries of the trajectory.
     */
    std::unordered_map<long, std::vector<Point<>>> grid_traj(point_list_t::const_iterator traj_b,
                                                            point_list_t::const_iterator traj_e,
                                                            double chord_l,
                                                            double& ux, double &uy, double& lx, double& ly) {


        auto last_pt = traj_b;

        if (last_pt == traj_e) {
            return {};
        }

        std::tie(lx, ly, ux, uy) = bounding_box(traj_b, traj_e);

        long g_size = static_cast<long>((ux - lx) / chord_l) + 1;
        std::unordered_map<long, std::vector<Point<>>> traj_points;
        if (traj_e - traj_b == 1) {
            //If this is a single point then we just return the point.
            traj_points.emplace(0, std::vector<Point<>>(traj_b, traj_e));
            return traj_points;
        }

        long location = index((*last_pt)(0), (*last_pt)(1), lx, ly, chord_l, g_size);
        traj_points.emplace(location, std::initializer_list<Point<>>{*last_pt});

        for (auto curr_pt = last_pt + 1; curr_pt != traj_e; curr_pt++) {

            auto g_x = static_cast<int>(((*last_pt)(0) - lx) / chord_l);
            auto g_y = static_cast<int>(((*last_pt)(1) - ly) / chord_l);
            auto g_n_x = static_cast<int>(((*curr_pt)(0) - lx) / chord_l);
            auto g_n_y = static_cast<int>(((*curr_pt)(1) - ly) / chord_l);

            if (g_n_x < g_x) {
                std::swap(g_n_x, g_x);
            }
            if (g_n_y < g_y) {
                std::swap(g_n_y, g_y);
            }

            for (int i = g_x; i <= g_n_x; i++) {
                double x_val = i * chord_l + lx;
                double y_val = x_to_y((*last_pt)(0), (*last_pt)(1), (*curr_pt)(0), (*curr_pt)(1), x_val);
                if (!inside_box(x_val, y_val, *last_pt, *curr_pt)) {
                    continue;
                }
                int j = static_cast<int>((y_val - ly) / chord_l);
                auto element = traj_points.find(i + g_size * j);
                if (element == traj_points.end()) {
                    traj_points.emplace(i + g_size * j, std::initializer_list<Point<>>{Point<>(x_val, y_val, 1.0)});
                } else {
                    element->second.emplace_back(x_val, y_val, 1.0);

                }
            }
            for (int j = g_y; j <= g_n_y; j++) {
                double y_val = j * chord_l + ly;
                double x_val = y_to_x((*last_pt)(0), (*last_pt)(1), (*curr_pt)(0), (*curr_pt)(1), y_val);
                if (!inside_box(x_val, y_val, *last_pt, *curr_pt)) {
                    continue;
                }
                int i = static_cast<int>((x_val - lx) / chord_l);
                auto element = traj_points.find(i + g_size * j);
                if (element == traj_points.end()) {
                    traj_points.emplace(i + g_size * j, std::initializer_list<Point<>>{Point<>(x_val, y_val, 1.0)});
                } else {
                    element->second.emplace_back(x_val, y_val, 1.0);
                }
            }

            location = index((*curr_pt)(0), (*curr_pt)(1), lx, ly, chord_l, g_size);
            auto element = traj_points.find(location);
            if (element == traj_points.end()) {
                traj_points.emplace(location, std::initializer_list<Point<>>{*curr_pt});
            } else {
                element->second.push_back(*curr_pt);
            }

            last_pt = curr_pt;
        }
        return traj_points;
    }


    /*
     * This is a wrapper for grid_traj that ensures all points are in their own cells.
     */
    point_list_t  grid_traj(point_list_t const& traj, double grid_resoluation) {
        point_list_t pts;
        double ly, lx, ux, uy;
        for (auto& elements : grid_traj(traj.begin(), traj.end(), grid_resoluation, ux, uy, lx, ly)) {
            pts.insert(pts.end(), elements.second.begin(), elements.second.end());
        }
        remove_duplicates(pts);
        return pts;
    }

    /*
    * Takes a trajectory and grids it so that each grid contains a single point of a trajectory that crosses it.
     * Chooses the cell to be in the center of each grid cell.
    */
    point_list_t  approx_traj_grid(point_list_t const& trajectory, double grid_resolution) {


        auto traj_b = trajectory.begin(), traj_e = trajectory.end();
        auto last_pt = traj_b;

        if (last_pt == traj_e) {
            return {};
        }

        double lx, ly, ux, uy;
        std::tie(lx, ly, ux, uy) = bounding_box(traj_b, traj_e);
        long g_size = static_cast<long>((ux - lx) / grid_resolution) + 1;
        std::unordered_set<long> traj_labels;

        long location = index((*last_pt)(0), (*last_pt)(1), lx, ly, grid_resolution, g_size);
        traj_labels.emplace(location);

        for (auto curr_pt = last_pt + 1; curr_pt != traj_e; curr_pt++) {

            auto g_x = static_cast<int>(((*last_pt)(0) - lx) / grid_resolution);
            auto g_y = static_cast<int>(((*last_pt)(1) - ly) / grid_resolution);
            auto g_n_x = static_cast<int>(((*curr_pt)(0) - lx) / grid_resolution);
            auto g_n_y = static_cast<int>(((*curr_pt)(1) - ly) / grid_resolution);

            if (g_n_x < g_x) {
                std::swap(g_n_x, g_x);
            }
            if (g_n_y < g_y) {
                std::swap(g_n_y, g_y);
            }

            for (int i = g_x + 1; i <= g_n_x; i++) {
                double x_val = i * grid_resolution + lx;
                double y_val = x_to_y((*last_pt)(0), (*last_pt)(1), (*curr_pt)(0), (*curr_pt)(1), x_val);
                if (!inside_box(x_val, y_val, *last_pt, *curr_pt)) {
                    continue;
                }
                int j = static_cast<int>((y_val - ly) / grid_resolution);
                auto element = traj_labels.find(i + g_size * j);
                if (element == traj_labels.end()) {
                    traj_labels.emplace(i + g_size * j);
                }
            }
            for (int j = g_y + 1; j <= g_n_y; j++) {
                double y_val = j * grid_resolution + ly;
                double x_val = y_to_x((*last_pt)(0), (*last_pt)(1), (*curr_pt)(0), (*curr_pt)(1), y_val);
                if (!inside_box(x_val, y_val, *last_pt, *curr_pt)) {
                    continue;
                }
                int i = static_cast<int>((x_val - lx) / grid_resolution);
                auto element = traj_labels.find(i + g_size * j);
                if (element == traj_labels.end()) {
                    traj_labels.emplace(i + g_size * j);
                }
            }

            location = index((*curr_pt)(0), (*curr_pt)(1), lx, ly, grid_resolution, g_size);
            auto element = traj_labels.find(location);
            if (element == traj_labels.end()) {
                traj_labels.emplace(location);
            }
            last_pt = curr_pt;
        }

        //Now convert the set into a set of points
        point_list_t simplified_traj;
        for (auto label : traj_labels) {
            long i = label % g_size;
            long j = label / g_size;
            //std::cout << i << " " << j << std::endl;
            double x_val_low = i * grid_resolution + lx;
            double y_val_low = j * grid_resolution + ly;
            double x_val_up = (i + 1) * grid_resolution + lx;
            double y_val_up = (j + 1) * grid_resolution + ly;
            simplified_traj.emplace_back((x_val_low + x_val_up) / 2.0, (y_val_up + y_val_low) / 2.0, 1.0);
        }
        return simplified_traj;
    }

    std::unordered_map<long, std::vector<Point<>>>
    approximate_traj_cells(point_list_t::const_iterator traj_b,
                            point_list_t::const_iterator traj_e, double chord_l, double eps) {
        double ux, uy, lx, ly;
        auto cells = grid_traj(traj_b, traj_e, chord_l, ux, uy, lx, ly);
        for (auto b = cells.begin(); b != cells.end(); b++) {
            auto approx = eps_core_set(eps, [&](Vec2 const& direction) {
                auto pt_max = std::max_element(b->second.begin(), b->second.end(),
                                               [&] (Point<> const& p1, Point<> const& p2){
                                                   double mag_1 = p1(0) * direction[0] + p1(1) * direction[1];
                                                   double mag_2 = p2(0) * direction[0] + p2(1) * direction[1];
                                                   return mag_1 < mag_2;
                                               });
                return Vec2{(*pt_max)(0), (*pt_max)(1)};
            });

            b->second.clear();
            for (auto approx_b = approx.begin(); approx_b != approx.end(); approx_b++) {
                b->second.emplace_back((*approx_b)[0], (*approx_b)[1], 1.0);
            }
        }
        return cells;
    }

    void approx_traj(point_list_t::const_iterator traj_b, point_list_t::const_iterator traj_e,
                        double chord_l, double eps, point_list_t& output) {
        for (auto& elements : approximate_traj_cells(traj_b, traj_e, chord_l, eps)) {
            output.insert(output.end(), elements.second.begin(), elements.second.end());
        }
    }


    void approx_traj_labels(point_list_t::const_iterator traj_b, point_list_t::const_iterator traj_e,
                     double chord_l, double eps, size_t label, double weight, lpoint_list_t& output) {
        for (auto& elements : approximate_traj_cells(traj_b, traj_e, chord_l, eps)) {
            for (auto& pt : elements.second) {
                output.emplace_back(label, weight, pt[0], pt[1], pt[2]);
            }
        }
        remove_duplicates(output);
    }

    point_list_t approx_traj_kernel_grid(point_list_t const& trajectory_pts, double chord_l, double eps) {

        point_list_t output;
        for (auto& elements : approximate_traj_cells(trajectory_pts.begin(), trajectory_pts.end(), chord_l, eps)) {
            for (auto& pt : elements.second) {
                //Create a point with the desired weight and label.
                output.emplace_back(pt[0], pt[1], pt[2]);
            }
        }
        remove_duplicates(output);
        return output;
    }

    point3_list_t kernel3d(point3_list_t const& pts, double eps) {

        auto glpts = new glReal[3 * pts.size()];

        for (size_t i = 0; i < pts.size(); i++) {
            glpts[i * 3    ] = pts[i](0);
            glpts[i * 3 + 1] = pts[i](1);
            glpts[i * 3 + 2] = pts[i](2);
        }
        glPointSet kernel_set;
        assert(pts.size() <= static_cast<size_t>(std::numeric_limits<int>::max()));

        kernel_set.init(3, static_cast<int>(pts.size()), glpts);


        glPointSet* p = kernel_set.getCoreSet3(static_cast<int>(2 / eps), eps / 2); // compute robust coreset

        point3_list_t core_set;
        // the reason why this starts at 1 is because the first point is always 0,0 for some reason.
        for (int i = 1; i < p->getNum(); i++) {
            glReal pt[3];
            p->getPoint(i, pt);
            core_set.emplace_back(pt[0], pt[1], pt[2], 1.0);
        }
        p->dump();
        delete[] glpts;
        kernel_set.dump();
        free(p);
        return core_set;
    }


    inline pt3_t lift_pt(pt2_t const&  pt) {
        double x = pt(0);
        double y = pt(1);
        return pt3_t(x, y, x * x + y * y, 1.0);
    }

    point_list_t lifting_coreset(point_list_t const& segments, double eps) {
        point_list_t pts = grid_traj(segments, eps / 2.0);
        point3_list_t lifted_set(pts.size(), pt3_t());
        std::transform(pts.begin(), pts.end(), lifted_set.begin(), lift_pt);
        auto pt3_set = kernel3d(lifted_set, eps);
        point_list_t  output_pts;
        for (auto& pt : pt3_set) {
            output_pts.emplace_back(pt(0), pt(1), 1.0);
        }
        return output_pts;
    }


    // This function does a Ramer-Douglas-Peucker Compression over a Trajectory
    // For more info, pls refer to:
    // https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm
    void rdp_compress(point_list_t::const_iterator begin, point_list_t::const_iterator end,
                     point_list_t& compressed, double eps) {
        if (end - begin <= 1) {
            return;
        }
        auto max_it = begin;
        double max = 0.0;
        for (auto it = begin + 1; it != end - 1; ++it) {
            auto tmp_dis = it->square_dist(*begin, *(end - 1));
            if (max < tmp_dis) {
                max = tmp_dis;
                max_it = it;
            }
        }
        if (max > eps * eps) {
            rdp_compress(begin, max_it + 1, compressed, eps);
            compressed.pop_back();
            rdp_compress(max_it, end, compressed, eps);
        } else {
            compressed.push_back(*begin);
            compressed.push_back(*(end - 1));
        }
    }

    point_list_t dp_compress(const point_list_t& trajectory, double eps) {
        point_list_t simplified_traj;
        rdp_compress(trajectory.begin(),  trajectory.end(), simplified_traj, eps);
        return simplified_traj;
    }

    double total_length(const point_list_t& pts) {
        return trajectory_t(pts).get_length();
    }

    point_list_t interval_sample(const trajectory_set_t& trajectories, std::vector<double>& indices, bool take_endpoints, bool take_single) {
        std::sort(indices.begin(), indices.end());

        double scaling_fact = std::accumulate(trajectories.begin(), trajectories.end(), 0.0,
                [&](const double& cum_weight, const trajectory_t& traj) {
            return cum_weight + traj.get_length() * traj.get_weight();
        });


        point_list_t sample_pts;
        auto idx = indices.begin();
        double curr_length = 0;
        for (auto& traj : trajectories) {
            if (traj.empty()) {
                continue;
            }
            auto last_pt = *(traj.begin());
            bool one_taken = take_endpoints;
            if (take_endpoints) {
                sample_pts.push_back(last_pt);
            }

            for (auto traj_b = traj.begin() + 1; traj_b != traj.end(); traj_b++) {

                double seg_length = (last_pt.dist(*traj_b) * traj.get_weight()) / scaling_fact;

                while (idx != indices.end() && seg_length + curr_length > *idx) {
                    double inside_length = (*idx - curr_length) / traj.get_weight();
                    double alpha = inside_length * scaling_fact / last_pt.dist(*traj_b);
                    alpha = std::min(1.0, std::max(alpha, 0.0));
                    // Create a new point on this line segment scaled between the two.
                    one_taken = true;
                    sample_pts.emplace_back(last_pt.on_segment(*traj_b, alpha));
                    idx++;
                }
                if (take_endpoints) {
                    sample_pts.push_back(*traj_b);
                }
                //Increment last_pt to the current point and update the total length.
                last_pt = *traj_b;
                curr_length += seg_length;
            }
            if (!one_taken && take_single) {
                sample_pts.push_back(last_pt);
            }
        }
        return sample_pts;
    }

    point_list_t interval_sample_single(const point_list_t& traj, std::vector<double>& indices, bool take_endpoints, bool take_single) {
        std::sort(indices.begin(), indices.end());

        double scaling_fact = trajectory_t(traj).get_length();

        point_list_t sample_pts;
        auto idx = indices.begin();
        double curr_length = 0;
        if (traj.empty()) {
            return sample_pts;
        }
        auto last_pt = *(traj.begin());
        bool one_taken = take_endpoints;
        if (take_endpoints) {
            sample_pts.push_back(last_pt);
        }

        for (auto traj_b = traj.begin() + 1; traj_b != traj.end(); traj_b++) {

            double seg_length = last_pt.dist(*traj_b) / scaling_fact;

            while (idx != indices.end() && seg_length + curr_length > *idx) {
                double inside_length = (*idx - curr_length);
                double alpha = inside_length * scaling_fact / last_pt.dist(*traj_b);
                alpha = std::min(1.0, std::max(alpha, 0.0));
                // Create a new point on this line segment scaled between the two.
                one_taken = true;
                sample_pts.emplace_back(last_pt.on_segment(*traj_b, alpha));
                idx++;
            }
            if (take_endpoints) {
                sample_pts.push_back(*traj_b);
            }
            //Increment last_pt to the current point and update the total length.
            last_pt = *traj_b;
            curr_length += seg_length;
        }
        if (!one_taken && take_single) {
            sample_pts.push_back(last_pt);
        }
        return sample_pts;
    }

    point_list_t uniform_sample(const trajectory_set_t& trajectories, size_t s, bool take_endpoints) {
        std::random_device rd;
        std::mt19937   generator(rd());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        std::vector<double> indices;
        for (size_t i = 0; i < s; i++) {
            indices.push_back(distribution(generator));
        }
        return interval_sample(trajectories, indices, take_endpoints, false);
    }

    point_list_t even_sample(const trajectory_set_t& trajectories, size_t s, bool take_endpoints) {
        /*
         * sample s points at 1/s distance appart.
         */
        std::vector<double> indices;
        for (size_t i = 0; i < s; i++) {
            indices.emplace_back( (i + .5) / s);
        }
        return interval_sample(trajectories, indices, take_endpoints, false);
    }

    point_list_t block_sample(const trajectory_set_t& trajectories, size_t s, bool take_endpoints) {
        /*
         * Chooses a single point uniformly at random in every alpha block of length L.
         */
        std::random_device rd;
        std::mt19937  generator(rd());
        std::uniform_real_distribution<double> distribution(0.0, 1.0 / s);
        std::vector<double> indices;
        for (size_t i = 0; i < s; i++) {
            indices.push_back(i / static_cast<double>(s) + distribution(generator));
        }
        return interval_sample(trajectories, indices, take_endpoints, false);
    }


    point_list_t uniform_sample_error(const point_list_t& traj, double eps, bool take_endpoints) {
        double wl = total_length(traj);
        assert(wl >= 0);
        auto s = static_cast<size_t>(std::lround(wl / eps));
        std::random_device rd;
        std::mt19937   generator(rd());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        std::vector<double> indices;
        for (size_t i = 0; i < s; i++) {
            indices.push_back(distribution(generator));
        }
        return interval_sample_single(traj, indices, take_endpoints, true);
    }

    point_list_t even_sample_error(const point_list_t& traj, double eps, bool take_endpoints){
        double wl = total_length(traj);
        assert(wl >= 0);
        auto s = static_cast<size_t>(std::lround(wl / eps));
        std::vector<double> indices;
        for (size_t i = 0; i < s; i++) {
            indices.emplace_back( (i + .5) / s);
        }
        return interval_sample_single(traj, indices, take_endpoints, true);
    }

    point_list_t block_sample_error(const point_list_t& traj, double eps, bool take_endpoints) {
        double wl = total_length(traj);
        assert(wl >= 0);
        auto s = static_cast<size_t>(std::lround(wl / eps));
        std::random_device rd;
        std::mt19937  generator(rd());
        std::uniform_real_distribution<double> distribution(0.0, 1.0 / s);
        std::vector<double> indices;
        for (size_t i = 0; i < s; i++) {
            indices.push_back(i / static_cast<double>(s) + distribution(generator));
        }
        return interval_sample_single(traj, indices, take_endpoints, true);
    }


    std::tuple<halfspace2_t, double> error_halfplane_coreset(const trajectory_t& trajectory, const point_list_t& pts) {
        halfspace2_t max_plane;
        double eps = 0;
        for (size_t i = 0; i < pts.size() - 1; i++) {
            for (size_t j = i + 1; j < pts.size(); j++) {
                halfspace2_t curr_plane(pts[i], pts[j]);

                auto norm_f = [&curr_plane](pt2_t const& p1, pt2_t const& p2) {
                    return p1.pdot(curr_plane.get_coords()) < p2.pdot(curr_plane.get_coords());
                };
                auto [min_el, max_el] = std::minmax_element(trajectory.begin(), trajectory.end(), norm_f);
                auto [min_aprx, max_aprx] = std::minmax_element(pts.begin(), pts.end(), norm_f);

                double eps_1 = std::abs(curr_plane.get_coords().pdot(*min_el) - curr_plane.get_coords().pdot(*min_aprx));
                double eps_2 = std::abs(curr_plane.get_coords().pdot(*max_el) - curr_plane.get_coords().pdot(*max_aprx));
                if (eps < eps_1 || eps < eps_2) {
                    eps = std::max(eps_1, eps_2);
                    max_plane = curr_plane;
                }
            }
        }
        return std::make_tuple(max_plane, eps);
    }

//    std::tuple<Disk, double> error_disk_coreset(const trajectory_t& trajectory,
//                double min_radius,
//                double max_radius,
//                const point_list_t& pts) {
//        Disk max_disk;
//        if (3 > pts.size()){
//            return std::make_tuple(max_disk, std::numeric_limits<double>::infinity());
//        }
//        double eps = 0;
//        for (size_t i = 0; i < pts.size() - 2; i++) {
//            for (size_t j = i + 1; j < pts.size() - 1; j++) {
//                for (size_t k = j + 1; k < pts.size(); k++) {
//                    Disk curr_disk(pts[i], pts[j], pts[k]);
//                    if (curr_disk.getRadius() > max_radius || min_radius > curr_disk.getRadius()) {
//                        continue;
//                    }
//                    auto min_dist = trajectory.point_dist(curr_disk.getOrigin());
//                    double error = std::abs(curr_disk.getRadius() - min_dist);
//                    if (eps < error) {
//                        eps = error;
//                        max_disk = curr_disk;
//                    }
//                }
//            }
//        }
//        return std::make_tuple(max_disk, eps);
//    }



}