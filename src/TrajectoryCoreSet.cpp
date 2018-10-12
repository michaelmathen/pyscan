//
// Created by mmath on 9/28/18.
//

#include <random>

#include "appext.h"
#include "Point.hpp"
#include "TrajectoryCoreSet.hpp"

namespace pyscan {




    /*
  * Compute the lower leftmost corner of a box containing these points.
  */
    std::tuple<double, double, double, double> bounding_box(point_it traj_b, point_it traj_e) {


        auto bounds_x = std::minmax_element(traj_b, traj_e,
                                            [&](Point<> const& p1, Point<> const& p2){
                                                return p1(0) < p2(0);
                                            });
        auto bounds_y = std::minmax_element(traj_b, traj_e, [&](Point<> const& p1, Point<> const& p2){
            return p1(1) <p2(1);
        });

        return std::make_tuple((*bounds_x.first)(0), (*bounds_y.first)(0),
                               (*bounds_x.second)(1), (*bounds_y.second)(1));
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

            for (int i = g_x + 1; i <= g_n_x; i++) {
                double x_val = i * chord_l + lx;
                double y_val = x_to_y((*last_pt)(0), (*last_pt)(1), (*curr_pt)(0), (*curr_pt)(1), x_val);
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
                double x_val = y_to_x((*last_pt)(0), (*last_pt)(1), (*curr_pt)(0), (*curr_pt)(1), y_val);
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


    std::unordered_map<long, std::vector<Point<>>>
    approximate_traj_cells(point_it traj_b, point_it traj_e, double chord_l, double eps) {
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

    void approx_traj(point_it traj_b, point_it traj_e, double chord_l, double eps, point_list_t& output) {
        for (auto& elements : approximate_traj_cells(traj_b, traj_e, chord_l, eps)) {
            output.insert(output.end(), elements.second.begin(), elements.second.end());
        }
    }

    point_list_t approx_traj_labels(point_list_t const& trajectory_pts, double weight, double chord_l, double eps) {

        point_list_t output;
        for (auto& elements : approximate_traj_cells(traj_b, traj_e, chord_l, eps)) {
            for (auto& pt : elements.second) {
                //Create a point with the desired weight and label.
                output.emplace_back(label, weight, pt[0], pt[1], pt[2]);
            }
        }
        return output;
    }




    const int KERNEL_LEVELS = 3;
    point3_list_t kernel3d(point3_list_t const& pts, double eps) {

        glReal* glpts = new glReal[3 * pts.size()];

        for (size_t i = 0; i < pts.size(); i++) {
            glpts[i * 3    ] = pts[i](0);
            glpts[i * 3 + 1] = pts[i](1);
            glpts[i * 3 + 2] = pts[i](2);
        }
        glPointSet kernel_set;
        kernel_set.init(3, pts.size(), glpts);


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

    point_list_t lifting_coreset(point_list_t const& pts, double eps) {
        point3_list_t lifted_set(pts.size(), pt3_t());
        std::transform(pts.begin(), pts.end(), lifted_set.begin(), lift_pt);
        auto pt3_set = kernel3d(lifted_set, eps);
        point_list_t  output_pts;
        for (auto& pt : pt3_set) {
            output_pts.emplace_back(pt(0), pt(1), 1.0);
        }
        return output_pts;
    }

}