//
// Created by mmath on 9/28/18.
//

#include <random>

#include "Point.hpp"

#include "gdiam.hpp"


namespace pyscan {


    point3_list_t hexagonal_packing(size_t pt_num) {

        double epsilon = 1e-5;

        //TODO generate two random points.
        point3_list_t sphere_pts;

        for (size_t i = 2; i < pt_num; i++) {
            while (true) {

                //Apply updates on all points
                for (auto& pt : sphere_pts) {

                    
                }
            }
        }
    }

    auto approx_mvbb(point_list_t const& pts, double eps) {

        auto all_pts = new gdiam_point[pts.size()];
        auto gpoint_mem = new gdiam_real[3 * pts.size()];

        for (size_t i = 0; i < pts.size(); i++) {
            gpoint_mem[3 * i] = pts[i](0);
            gpoint_mem[3 * i + 1] = pts[i](1);
            gpoint_mem[3 * i + 2] = pts[i](2);
            all_pts[i] = &gpoint_mem[3 * i];
        }
        auto bbox = gdiam_approx_mvbb(all_pts, pts.size(), eps);

    }


}