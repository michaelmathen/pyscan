//
// Created by mmath on 9/4/18.
//
#include <tuple>
#include <unordered_map>

#include "FunctionApprox.hpp"
#include "TrajectoryScan.hpp"


namespace pyscan {

    /*
     * Compute the lower leftmost corner of a box containing these points.
     */
    std::tuple<double, double, double, double> bounding_box(point_it traj_b, point_it traj_e) {


        auto bounds_x = std::minmax_element(traj_b, traj_e,
                [&](Point<> const& p1, Point<> const& p2){
                    return getX(p1) < getX(p2);
        });
        auto bounds_y = std::minmax_element(traj_b, traj_e, [&](Point<> const& p1, Point<> const& p2){
            return getY(p1) < getY(p2);
        });

        return std::make_tuple(getX(*bounds_x.first), getY(*bounds_y.first),
                getX(*bounds_x.second), getY(*bounds_y.second));
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
    /*
     * Takes a trajectory and grids it so that each grid contains points that cross it..
     */
    std::unordered_map<long, std::vector<Point<>>>
            grid_traj(point_it traj_b, point_it traj_e, double chord_l) {


        auto last_pt = traj_b;

        if (last_pt == traj_e) {
            return {};
        }

        double lx, ly, ux, uy;
        std::tie(lx, ly, ux, uy) = bounding_box(traj_b, traj_e);

        long g_size = static_cast<long>((ux - lx) / chord_l) + 1;
        std::unordered_map<long, std::vector<Point<>>> traj_points;

        long location = index(getX(*last_pt), getY(*last_pt), lx, ly, chord_l, g_size);
        traj_points.emplace(location, std::initializer_list<Point<>>{*last_pt});

        for (auto curr_pt = last_pt + 1; curr_pt != traj_e; curr_pt++) {

            auto g_x = static_cast<int>((getX(*last_pt) - lx) / chord_l);
            auto g_y = static_cast<int>((getY(*last_pt) - ly) / chord_l);
            auto g_n_x = static_cast<int>((getX(*curr_pt) - lx) / chord_l);
            auto g_n_y = static_cast<int>((getY(*curr_pt) - ly) / chord_l);

            if (g_n_x < g_x) {
                std::swap(g_n_x, g_x);
            }
            if (g_n_y < g_y) {
                std::swap(g_n_y, g_y);
            }

            for (int i = g_x + 1; i <= g_n_x; i++) {
                double x_val = i * chord_l + lx;
                double y_val = x_to_y(getX(*last_pt), getY(*last_pt), getX(*curr_pt), getY(*curr_pt), x_val);
                int j = static_cast<int>((y_val - ly) / chord_l);
                auto element = traj_points.find(i + g_size * j);
                if (element == traj_points.end()) {
                    traj_points.emplace(i + g_size * j, std::initializer_list<Point<>>{Point<>(x_val, y_val, 1.0)});
                } else {
                    element->second.emplace_back(x_val, y_val, 1.0);

                }
            }
            for (int j = g_y + 1; j <= g_n_y; j++) {
                double y_val = j * chord_l + ly;
                double x_val = y_to_x(getX(*last_pt), getY(*last_pt), getX(*curr_pt), getY(*curr_pt), y_val);
                int i = static_cast<int>((x_val - lx) / chord_l);
                auto element = traj_points.find(i + g_size * j);
                if (element == traj_points.end()) {
                    traj_points.emplace(i + g_size * j, std::initializer_list<Point<>>{Point<>(x_val, y_val, 1.0)});
                } else {
                    element->second.emplace_back(x_val, y_val, 1.0);
                }
            }

            location = index(getX(*curr_pt), getY(*curr_pt), lx, ly, chord_l, g_size);
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


    std::unordered_map<long, std::vector<Point<>>>
            approximate_traj_cells(point_it traj_b, point_it traj_e, double chord_l, double eps) {
        auto cells = grid_traj(traj_b, traj_e, chord_l);
        for (auto b = cells.begin(); b != cells.end(); b++) {
            auto approx = eps_core_set(eps, [&](Vec2 const& direction) {
                auto pt_max = std::max_element(b->second.begin(), b->second.end(),
                        [&] (Point<> const& p1, Point<> const& p2){
                    double mag_1 = getX(p1) * direction[0] + getY(p1) * direction[1];
                    double mag_2 = getX(p2) * direction[0] + getY(p2) * direction[1];
                    return mag_1 < mag_2;
                });
                return Vec2{getX(*pt_max), getY(*pt_max)};
            });

            b->second.clear();
            for (auto approx_b = approx.begin(); approx_b != approx.end(); approx_b++) {
                b->second.emplace_back((*approx_b)[0], (*approx_b)[1], 1.0);
            }
        }
        return cells;
    }

    void approx_traj(point_it traj_b, point_it traj_e, double chord_l, double eps, point_list& output) {
        for (auto& elements : approximate_traj_cells(traj_b, traj_e, chord_l, eps)) {
            output.insert(output.end(), elements.second.begin(), elements.second.end());
        }
    }

    void approx_traj_labels(point_it traj_b, point_it traj_e, size_t label,
            double weight, double chord_l, double eps, lpoint_list& output) {
        for (auto& elements : approximate_traj_cells(traj_b, traj_e, chord_l, eps)) {
            output.reserve(output.size() + elements.second.size());
            for (auto& pt : elements.second) {
                //Create a point with the desired weight and label.
                output.emplace_back(label, weight, pt[0], pt[1], pt[2]);
            }
        }
    }



    std::vector<Point<>> core_set_3d_traj(point_list traj, double alpha) {

        auto approx = eps_core_set3(alpha, [&](Vec3 const& direction) {
            auto pt_max = std::max_element(traj.begin(), traj.end(), [&] (Point<> const& p1, Point<> const& p2){
                double z_1 = getX(p1) * getX(p1) + getY(p1) * getY(p1);
                double z_2 = getX(p2) * getX(p2) + getY(p2) * getY(p2);

                double mag_1 = getX(p1) * direction[0] + getY(p1) * direction[1] + z_1 * direction[2];
                double mag_2 = getX(p2) * direction[0] + getY(p2) * direction[1] + z_2 * direction[2];
                return mag_1 < mag_2;
            });
            return Vec3{getX(*pt_max), getY(*pt_max), getX(*pt_max) * getX(*pt_max) + getY(*pt_max) * getY(*pt_max)};
        });

        std::vector<Point<>> output;
        for (auto& pt : approx) {
            // Do a flat projection back down.
            output.emplace_back(pt[0], pt[1], 1.0);
        }
        return output;
    }

    void core_set_3d_traj_internal(point_it traj_b, point_it traj_e, size_t label, double weight, double alpha,
                          lpoint_list& output) {

        auto approx = eps_core_set3(alpha, [&](Vec3 const& direction) {
            auto pt_max = std::max_element(traj_b, traj_e, [&] (Point<> const& p1, Point<> const& p2){
                double z_1 = getX(p1) * getX(p1) + getY(p1) * getY(p1);
                double z_2 = getX(p2) * getX(p2) + getY(p2) * getY(p2);

                double mag_1 = getX(p1) * direction[0] + getY(p1) * direction[1] + z_1 * direction[2];
                double mag_2 = getX(p2) * direction[0] + getY(p2) * direction[1] + z_2 * direction[2];
                return mag_1 < mag_2;
            });
            return Vec3{getX(*pt_max), getY(*pt_max), getX(*pt_max) * getX(*pt_max) + getY(*pt_max) * getY(*pt_max)};
        });

        for (auto& pt : approx) {
            output.emplace_back(label, weight, pt[0], pt[1], 1.0);
        }
    }
    std::tuple<Disk, double> traj_disk_scan(traj_set &net,
                                            wtraj_set &sampleM,
                                            wtraj_set &sampleB,
                                            double alpha,
                                            double min_r,
                                            std::function<double(double, double)> const &scan) {


        double chord_l = sqrt(4 * alpha * min_r - 2 * alpha * alpha);
        // Go through each trajectory and approximate it with a coreset of points

        // Go through the set of net points.
        point_list net_points;
        size_t offset = 0;
        for(auto b = net.offsets.begin(); b != net.offsets.end(); b++) {
            auto traj_b = net.traj_pts.begin() + offset;
            auto traj_e = net.traj_pts.begin() + *b;
            offset = *b;
            approx_traj(traj_b, traj_e, chord_l, alpha, net_points);
        }
        // go through the set of measured points
        lpoint_list sampleM_points;
        offset = 0;
        size_t label = 0;
        auto wb = sampleM.weights.begin();
        for(auto b = sampleM.offsets.begin(); b != sampleM.offsets.end(); b++) {
            auto traj_b = sampleM.traj_pts.begin() + offset;
            auto traj_e = sampleM.traj_pts.begin() + *b;
            offset = *b;
            approx_traj_labels(traj_b, traj_e, label, *wb, chord_l, alpha, sampleM_points);
            wb++;
            label++;//increment label
        }

        // go through the set of baseline points
        lpoint_list sampleB_points;
        offset = 0;
        label = 0;
        wb = sampleB.weights.begin();
        for(auto b = sampleB.offsets.begin(); b != sampleB.offsets.end(); b++) {
            auto traj_b = sampleB.traj_pts.begin() + offset;
            auto traj_e = sampleB.traj_pts.begin() + *b;
            offset = *b;
            approx_traj_labels(traj_b, traj_e, label, *wb, chord_l, alpha, sampleB_points);
            wb++;
            label++;//increment label
        }

        // Scan the resulting set of points using standard labeled disk scanning function.

        // return the max disk.
        auto grid_r = static_cast<uint32_t>(1.0 / min_r);
        std::cout << net_points.size() << " " << sampleM_points.size() << " " << sampleB_points.size() << std::endl;
        return disk_scan_scale(net_points, sampleM_points, sampleB_points, grid_r, scan);
    }

}