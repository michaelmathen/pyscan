//
// Created by mmath on 9/28/18.
//

#include <random>
#include <unordered_map>
#include <unordered_set>
#include <cassert>

#include "appext.h"
#include "Point.hpp"
#include "FunctionApprox.hpp"

#include "TrajectoryCoreSet.hpp"

namespace pyscan {



    template <typename T>
    void remove_duplicates(T& pts) {
        std::sort(pts.begin(), pts.end(), [](pt2_t const& p1, pt2_t const& p2){
            return p1(0) < p2(0);
        });

        auto end_it = std::unique(pts.begin(), pts.end(), [] (pt2_t const& p1, pt2_t const& p2) {
            return p1.approx_eq(p2);
        });
        pts.erase(end_it, pts.end());
    }


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
        return ux > x && x > lx &&  uy > y && y > ly;
    }
    /*
     * Takes a trajectory and grids it so that each grid contains points that cross it..
     */
    std::unordered_map<long, std::vector<Point<>>> grid_traj(point_list_t::const_iterator traj_b,
                                                            point_list_t::const_iterator traj_e,
                                                            double chord_l) {


        auto last_pt = traj_b;

        if (last_pt == traj_e) {
            return {};
        }

        double lx, ly, ux, uy;
        std::tie(lx, ly, ux, uy) = bounding_box(traj_b, traj_e);

        long g_size = static_cast<long>((ux - lx) / chord_l) + 1;
        std::unordered_map<long, std::vector<Point<>>> traj_points;

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


    point_list_t  grid_traj(point_list_t const& traj, double grid_resoluation) {
        point_list_t pts;

        for (auto& elements : grid_traj(traj.begin(), traj.end(), grid_resoluation)) {
            pts.insert(pts.end(), elements.second.begin(), elements.second.end());
        }
        remove_duplicates(pts);
        return pts;
    }

    /*
    * Takes a trajectory and grids it so that each grid contains a single point of a trajectory that crosses it.
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
        auto cells = grid_traj(traj_b, traj_e, chord_l);
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



    const int KERNEL_LEVELS = 3;
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
        return {x, y, x * x + y * y, 1.0};
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

}