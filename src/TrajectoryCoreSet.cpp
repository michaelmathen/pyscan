//
// Created by mmath on 9/28/18.
//

#include <random>

#include "appext.h"
#include "Point.hpp"
#include "TrajectoryCoreSet.hpp"

namespace pyscan {

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


        glPointSet* p = kernel_set.getRobustCoreset(static_cast<int>(2 / eps) , eps / 2, KERNEL_LEVELS); // compute robust coreset

        point3_list_t core_set;
        for (int i = 0; i < p->getNum(); i++) {
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