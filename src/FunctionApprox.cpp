//
// Created by mmath on 7/7/17.
//

#include <cmath>
#include <deque>
#include <algorithm>
#include <functional>
#include <tuple>
#include "Vecky.hpp"
#include "Utilities.hpp"

namespace pyscan {

    bool lineIntersection(VecD const& dir1, double d1, VecD const& dir2, double d2, VecD& v_ext) {
        /*
         * Need to test two cases.
         *
         */

        double mx1 = dir1[0] * d1 - dir2[0] * d2;
        double bx1 = -dir1[1] * d1 + dir2[1] * d2;
        double proj = (mx1 * dir1[0] + bx1 * dir1[1]);

        if (mx1 == 0 && bx1 == 0) {
            return false;
        } else if (std::signbit(d1) * proj >= 0) {
            v_ext[0] = mx1, v_ext[1] = bx1;
            return true;
        } else {
            v_ext[0] = -mx1, v_ext[1] = -bx1;
            return true;
        }
    }


    double approximateHull(double eps,
                           VecD const& cc, VecD const& cl,
                           std::function<double(VecD)> phi, //function to maximize
                           std::function<VecD(VecD)> lineMaxF) {


        auto avg = [&] (VecD const& v1, VecD const& v2) {
            VecD v_out = v1 + v2;
            return v_out * 1.0 / mag(v_out);
        };

        int ux, uy, lx, ly;
        struct Frame {
            VecD d_cc, d_cl, p_cc, p_cl;
            Frame(VecD const& vi, VecD const& vj, VecD const& cc, VecD const& cl) :
                    d_cc(vi), d_cl(vj), p_cc(cc), p_cl(cl) {}
        };
        double maxRValue = 0;

        std::deque<Frame> frameStack;
        // TODO double check debug to see if there is an issue with an infinite singularity here. Might need to change start.
        //This start needs to be fixed. Compute the mi, bi, and everything explicitly
        frameStack.push_back(Frame(cc, cl, lineMaxF(cc), lineMaxF(cl)));
        //cout << fixed;
        //cout << setprecision(9);
        while(!frameStack.empty()) {
            Frame lf = frameStack.front();
            frameStack.pop_front();
            double di = dot(lf.d_cc, lf.p_cc);
            double dj = dot(lf.d_cl, lf.p_cl);

            VecD p_ext;
            if (lineIntersection(lf.d_cc, di, lf.d_cl, dj, p_ext)) {
                continue;
            }
            double vw = phi(p_ext);
            double vi = phi(lf.p_cc);
            double vj = phi(lf.p_cl);

            if (std::max({vi, vj, maxRValue}) + eps  > vw) {
                // No value in this triangle is large enough to be interesting
                continue;
            } else {
                //This triangle is worth evaluating
                auto m_vec = avg(lf.d_cc, lf.d_cl);
                auto line_max = lineMaxF(m_vec);
                frameStack.push_back(Frame(m_vec, lf.d_cl,
                                           line_max, lf.p_cl));
                frameStack.push_back(Frame(lf.d_cc, m_vec,
                                           lf.p_cc, line_max));
            }
        }
        return maxRValue;
    }

    double approximateHull(double eps,
                           std::function<double(VecD)> phi, //function to maximize
                           std::function<VecD(VecD)> lineMaxF) {
        //This top line doesn't make sense to maximize for since we will always just return the largest region.
        //approximateHull(eps, VecD(0, 1), VecD(1, 0), f, linemaxF);
        //Likewise this line doesn't make any sense to optimize for since we will always try to find the null set.
        // approximateHull(eps, VecD(0, -1), VecD(-1, 0), f, linemaxF);

        //This finds the region with the largest measured amount, but smallest baseline amount in terms of the stat fun
        return std::max(approximateHull(eps, VecD(1, 0), VecD(0, -1), phi, lineMaxF),
                        //This finds the region with the largest baseline amount, but smallest measured amount in terms of the stat fun
                        approximateHull(eps, VecD(-1, 0), VecD(0, 1), phi, lineMaxF));
    }
  }
