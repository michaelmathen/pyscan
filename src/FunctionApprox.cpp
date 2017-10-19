//
// Created by mmath on 7/7/17.
//

#include <cmath>
#include <deque>
#include <algorithm>
#include <functional>
#include <tuple>
#include <iostream>

#include "Vecky.hpp"
#include "Utilities.hpp"

namespace pyscan {

    bool almostEqualRelative(double A, double B, double maxRelDiff) {
        // Calculate the difference.
        double diff = fabs(A - B);
        A = fabs(A);
        B = fabs(B);
        // Find the largest
        double largest = (B > A) ? B : A;
        return diff <= largest * maxRelDiff;
    }

    bool lineIntersection(VecD const& dir1, double d1, VecD const& dir2, double d2, VecD& v_ext) {
        /*
         * Need to test two cases.
         *
         */
        if (almostEqualRelative(dir1[0] * dir2[1], dir1[1] * dir2[0], 8 * std::numeric_limits<double>::epsilon())) {
            return false;
        }

        double denom = 1 / (dir1[0] * dir2[1] - dir1[1] * dir2[0]);

        v_ext[0] = (dir2[1] * d1 - dir1[1] * d2) * denom;
        v_ext[1] = (dir1[0] * d2 - dir2[0] * d1) * denom;
        return true;
    }

    bool inside(VecD pos, double alpha, double rho) {
      return (alpha <= pos[0]) && (pos[0] <= 1 - alpha) &&
        (rho <= pos[1]) && (pos[1] <= 1 - rho);
    }

    bool inRange(double a, double b) {
      return (b <= a) && (a <= 1 - b);
    }

    void lineIntersectionI(VecD const& dir1, double d1, VecD const& dir2, double d2, VecD& v_ext) {
      /*
        If there is numerical instability then we project to infinity.
      */
      if (!lineIntersection(dir1, d1, dir2, d2, v_ext)) {
        v_ext[0] = std::numeric_limits<double>::infinity();
        v_ext[1] = std::numeric_limits<double>::infinity();
      }
    }



    VecD projectToBoundary(VecD pos, VecD dir, double alpha, double rho) {
      /*
        Takes the pos and projects it to the rectangle defined by
        [rho, 1 rho] x [alpha, 1 - alpha]
      */
      if (inside(pos, alpha, rho))
        return pos;

      VecD rho_b, rho_t, alpha_l, alpha_r;
      double d1 = dot(pos, dir);
      lineIntersectionI(dir, d1, VecD(0.0, 1.0), rho, rho_b);
      lineIntersectionI(dir, d1, VecD(0.0, 1.0), 1 - rho, rho_t);
      lineIntersectionI(dir, d1, VecD(1.0, 0.0), alpha, alpha_l);
      lineIntersectionI(dir, d1, VecD(1.0, 0.0), 1 - alpha, alpha_r);
      if (pos[1] < rho && inRange(rho_b[0], alpha))
          return rho_b;
      else if (1 - rho < pos[1] && inRange(rho_t[0], alpha))
        return rho_t;
      else if (pos[0] < alpha && inRange(alpha_l[1], rho))
        return alpha_l;
      else if (1 - alpha < pos[0] && inRange(alpha_r[1], rho))
        return alpha_r;
      // lives in the corner and no good projections
      else if (pos[0] < alpha) {
        if (pos[1] < rho)
          return {alpha, rho};
        else
          return {alpha, 1 - rho};
      } else {
        if (pos[1] < rho)
          return {1 - alpha, rho};
        else
          return {1 - alpha, 1 - rho};
      }
    }

    double approximateHull(double eps,
                           VecD const& cc, VecD const& cl,
                           std::function<double(VecD)> phi, //function to maximize
                           std::function<VecD(VecD)> lineMaxExt) {


        auto lineMaxF = [&] (VecD v1) {
          auto pt = lineMaxExt(v1);
          //std::cout << "line_max = " << pt << std::endl;
          return pt; // projectToBoundary(pt, v1, alpha, rho);
        };

        auto avg = [&] (VecD const& v1, VecD const& v2) {
            VecD v_out = v1 + v2;
            VecD tmp;
            double norm = 1.0 / sqrt(v_out[0] * v_out[0] + v_out[1] * v_out[1]);
            tmp[0] = v_out[0] * norm;
            tmp[1] = v_out[1] * norm;
            return tmp;
        };

        int ux, uy, lx, ly;
        struct Frame {
            VecD d_cc, d_cl, p_cc, p_cl;
            Frame(VecD const& di, VecD const& dj, VecD const& cc, VecD const& cl) :
                    d_cc(di), d_cl(dj), p_cc(cc), p_cl(cl) {}
        };
        double maxRValue = 0;

        std::deque<Frame> frameStack;
        // TODO double check debug to see if there is an issue with an infinite singularity here. Might need to change start.
        //This start needs to be fixed. Compute the mi, bi, and everything explicitly
        frameStack.emplace_back(cc, cl, lineMaxF(cc), lineMaxF(cl));

//#ifndef NDEBUG
//        size_t iter_count = 0;
//#endif
        while(!frameStack.empty()) {
            Frame lf = frameStack.front();

//#ifndef NDEBUG
//            std::cout << "considering triangle" << std::endl;
//            std::cout << lf.d_cc << " " << lf.p_cc << " " << lf.d_cl << " " << lf.p_cl << std::endl;
//#endif
            frameStack.pop_front();
            double di = dot(lf.d_cc, lf.p_cc);
            double dj = dot(lf.d_cl, lf.p_cl);
            double vi = phi(lf.p_cc);
            double vj = phi(lf.p_cl);
            maxRValue = std::max({vi, vj, maxRValue});
            VecD p_ext;
            if (lineIntersection(lf.d_cc, di, lf.d_cl, dj, p_ext)) {
                double vw = phi(p_ext);
//#ifndef NDEBUG
//                std::cout << lf.p_cc << " " << vi << std::endl;
//                std::cout << p_ext << " " << vw << std::endl;
//                std::cout << lf.p_cl << " " << vj << std::endl;
//#endif
                if (maxRValue + eps <= vw) {
                    //This triangle is worth evaluating
//#ifndef NDEBUG
//                    std::cout << "evaluating triangle" << std::endl;
//#endif
                    VecD m_vec = avg(lf.d_cc, lf.d_cl);
                    auto line_max = lineMaxF(m_vec);

                    frameStack.emplace_back(lf.d_cc, m_vec, lf.p_cc, line_max);
                    frameStack.emplace_back(m_vec, lf.d_cl, line_max, lf.p_cl);
                }
            }
//#ifndef NDEBUG
//            iter_count += 1;
//            std::cout << std::endl;
//#endif
        }
//#ifndef NDEBUG
//        std::cout << "maxRValue = " << maxRValue << std::endl;
//        std::cout << iter_count << std::endl;
//#endif
        return maxRValue;
    }

    double approximateHull(double alpha, double rho, double eps,
                           std::function<double(VecD)> phi, //function to maximize
                           std::function<VecD(VecD)> lineMaxF) {
        //This top line doesn't make sense to maximize for since we will always just return the largest region.
        //approximateHull(eps, VecD(0, 1), VecD(1, 0), f, linemaxF);
        //Likewise this line doesn't make any sense to optimize for since we will always try to find the null set.
        // approximateHull(eps, VecD(0, -1), VecD(-1, 0), f, linemaxF);

        //This finds the region with the largest measured amount, but smallest baseline amount in terms of the stat fun
        return std::max(approximateHull(alpha, rho, eps, VecD(1, 0), VecD(0, -1), phi, lineMaxF),
                        //This finds the region with the largest baseline amount, but smallest measured amount in terms of the stat fun
                        approximateHull(alpha, rho, eps, VecD(-1, 0), VecD(0, 1), phi, lineMaxF));
    }
  }
