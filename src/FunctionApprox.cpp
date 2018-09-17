//
// Created by mmath on 7/7/17.
//

#include <cmath>
#include <deque>
#include <algorithm>
#include <functional>
#include <tuple>
#include <iostream>

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "Utilities.hpp"
#include "FunctionApprox.hpp"

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

    bool lineIntersection(Vec2 const& dir1, double d1, Vec2 const& dir2, double d2, Vec2& v_ext) {
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

    double det2(double a1, double a2, double b1, double b2) {
        return a1 * b2 - a2 * b1;
    }

    double det3(Vec3 const& dir1, Vec3 const& dir2, Vec3 const& dir3) {
        return dir1[0] * det2(dir2[1], dir2[2], dir3[1], dir3[2])
               - dir1[1] * det2(dir2[0], dir2[2], dir3[0], dir3[2])
               + dir1[2] * det2(dir2[0], dir2[1], dir3[0], dir3[1]);
    }

    bool lineIntersection(Vec3 const& dir1, double d1,
                          Vec3 const& dir2, double d2,
                          Vec3 const& dir3, double d3,
                          Vec3& v_ext) {
        /*
         * Need to test two cases.
         *
         */


        double dval = det3(dir1, dir2, dir3);
        if (fabs(dval) < 8 * std::numeric_limits<double>::epsilon()) {
            return false;
        }

        using namespace boost::numeric::ublas;

        matrix<double> A(3, 3);

        A(0, 0) = dir1[0], A(0, 1) = dir1[1], A(0, 2) = dir1[2],
        A(1, 0) = dir2[0], A(1, 1) = dir2[1], A(1, 2) = dir2[2],
        A(2, 0) = dir3[0], A(2, 1) = dir3[1], A(2, 2) = dir3[2];

        vector<double> x;
        x[0] = d1, x[1] = d2, x[2] = d3;
        permutation_matrix<size_t> P(3);
        lu_factorize(A, P);
        lu_substitute(A, P, x);

        v_ext[0] = x[0], v_ext[1] = x[1], v_ext[2] = x[2];
        return true;
    }


    bool inside(Vec2 const& pos, double alpha, double rho) {
      return (alpha <= pos[0]) && (pos[0] <= 1 - alpha) &&
        (rho <= pos[1]) && (pos[1] <= 1 - rho);
    }

    bool inRange(double a, double b) {
      return (b <= a) && (a <= 1 - b);
    }

    void lineIntersectionI(Vec2 const& dir1, double d1, Vec2 const& dir2, double d2, Vec2& v_ext) {
      /*
        If there is numerical instability then we project to infinity.
      */
      if (!lineIntersection(dir1, d1, dir2, d2, v_ext)) {
        v_ext[0] = std::numeric_limits<double>::infinity();
        v_ext[1] = std::numeric_limits<double>::infinity();
      }
    }



    Vec2 projectToBoundary(Vec2 const& pos, Vec2 const& dir, double alpha, double rho) {
      /*
        Takes the pos and projects it to the rectangle defined by
        [rho, 1 rho] x [alpha, 1 - alpha]
      */
      if (inside(pos, alpha, rho))
        return pos;

      Vec2 rho_b, rho_t, alpha_l, alpha_r;
      double d1 = dot(pos, dir);
      lineIntersectionI(dir, d1, Vec2{0.0, 1.0}, rho, rho_b);
      lineIntersectionI(dir, d1, Vec2{0.0, 1.0}, 1 - rho, rho_t);
      lineIntersectionI(dir, d1, Vec2{1.0, 0.0}, alpha, alpha_l);
      lineIntersectionI(dir, d1, Vec2{1.0, 0.0}, 1 - alpha, alpha_r);
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
                           Vec2 const& cc, Vec2 const& cl,
                           std::function<double(Vec2)> phi, //function to maximize
                           std::function<Vec2(Vec2)> lineMaxExt) {


        auto lineMaxF = [&] (Vec2 v1) {
          auto pt = lineMaxExt(v1);
          //std::cout << "line_max = " << pt << std::endl;
          return pt; // projectToBoundary(pt, v1, alpha, rho);
        };

        auto avg = [&] (Vec2 const& v1, Vec2 const& v2) {
            Vec2 v_out = v1 + v2;
            Vec2 tmp;
            double norm = 1.0 / sqrt(v_out[0] * v_out[0] + v_out[1] * v_out[1]);
            tmp[0] = v_out[0] * norm;
            tmp[1] = v_out[1] * norm;
            return tmp;
        };

        int ux, uy, lx, ly;
        struct Frame {
            Vec2 d_cc, d_cl, p_cc, p_cl;
            Frame(Vec2 const& di, Vec2 const& dj, Vec2 const& cc, Vec2 const& cl) :
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
            Vec2 p_ext;

            if (lineIntersection(lf.d_cc, di, lf.d_cl, dj, p_ext)) {
                double vw = phi(p_ext);
//#ifndef NDEBUG
//                std::cout << lf.p_cc << " " << vi << std::endl;
//                std::cout << p_ext << " " << vw << std::endl;
//                std::cout << lf.p_cl << " " << vj << std::endl;
//#endif
                //dist(lf)
                if (vw - maxRValue > eps) {

                    //This triangle is worth evaluating
//#ifndef NDEBUG
//                    std::cout << "evaluating triangle" << std::endl;
//#endif
                    Vec2 m_vec = avg(lf.d_cc, lf.d_cl);
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

    double approximateHull(double eps,
                           std::function<double(Vec2)> phi, //function to maximize
                           std::function<Vec2(Vec2)> lineMaxF) {
        return std::max(approximateHull(eps, Vec2{1, 0}, Vec2{0, -1}, phi, lineMaxF),
                        approximateHull(eps, Vec2{-1, 0}, Vec2{0, 1}, phi, lineMaxF));
    }

    double height(Vec2 const& p1, Vec2 const& p2, Vec2 const& pt) {
        double a = sqrt((p1[0] - p2[0]) *(p1[0] - p2[0]) + (p1[1] - p2[1]) *(p1[1] - p2[1]));
        double b = sqrt((p1[0] - pt[0]) *(p1[0] - pt[0]) + (p1[1] - pt[1]) *(p1[1] - pt[1]));
        double c = sqrt((p2[0] - pt[0]) *(p2[0] - pt[0]) + (p2[1] - pt[1]) *(p2[1] - pt[1]));
        double s = (a + b + c)/ 2;
        return 2 * sqrt(s * (s - a) * (s - b) * (s - c)) / a;
    }

    std::vector<Vec2> eps_core_set(double eps,
                                   Vec2 const& cc, Vec2 const& cl,
                                   std::function<Vec2(Vec2)> lineMaxF) {

            auto avg = [&] (Vec2 const& v1, Vec2 const& v2) {
                Vec2 v_out = v1 + v2;
                Vec2 tmp;
                double norm = 1.0 / sqrt(v_out[0] * v_out[0] + v_out[1] * v_out[1]);
                tmp[0] = v_out[0] * norm;
                tmp[1] = v_out[1] * norm;
                return tmp;
            };

            int ux, uy, lx, ly;
            struct Frame {
                Vec2 d_cc, d_cl, p_cc, p_cl;
                Frame(Vec2 const& di, Vec2 const& dj, Vec2 const& cc, Vec2 const& cl) :
                        d_cc(di), d_cl(dj), p_cc(cc), p_cl(cl) {}
            };
            auto pcc = lineMaxF(cc), pcl = lineMaxF(cl);
            std::vector<Vec2> pts{ pcc, pcl };
            std::deque<Frame> frameStack;
            frameStack.emplace_back(cc, cl, pcc, pcl);
            while(!frameStack.empty()) {
                Frame lf = frameStack.front();
                frameStack.pop_front();
                double di = dot(lf.d_cc, lf.p_cc);
                double dj = dot(lf.d_cl, lf.p_cl);

                Vec2 p_ext;

                //TODO fix this.
//                if (lineIntersection(lf.d_cc, di, lf.d_cl, dj, p_ext)) {
//                    if (height(lf.p_cc, lf.p_cl, p_ext) > eps) {
//
//                        Vec2 m_vec = avg(lf.d_cc, lf.d_cl);
//                        auto line_max = lineMaxF(m_vec);
//                        pts.push_back(line_max);
//                        frameStack.emplace_back(lf.d_cc, m_vec, lf.p_cc, line_max);
//                        frameStack.emplace_back(m_vec, lf.d_cl, line_max, lf.p_cl);
//                    }
//                }
            }
            return pts;
    }

    std::vector<Vec2> eps_core_set(double eps,
                        std::function<Vec2(Vec2)> lineMaxF) {
        auto core_set1 = eps_core_set(eps, Vec2{1, 0}, Vec2{0, -1},  lineMaxF)
            ,core_set2 = eps_core_set(eps, Vec2{0, 1}, Vec2{1, 0}, lineMaxF)
            ,core_set3 = eps_core_set(eps, Vec2{1, 0}, Vec2{-1, 0}, lineMaxF)
            ,core_set4 = eps_core_set(eps, Vec2{-1, 0}, Vec2{0, -1}, lineMaxF);
        core_set1.insert(core_set1.end(), core_set2.begin(), core_set2.end());
        core_set1.insert(core_set1.end(), core_set3.begin(), core_set3.end());
        core_set1.insert(core_set1.end(), core_set4.begin(), core_set4.end());
        return core_set1;
    }

    double height(Vec3 const& p1, Vec3 const& p2, Vec3 const& p3, Vec3 const& p_ext) {
        
    }

    std::vector<Vec3> eps_core_set3(double eps,
                                   Vec3 const& cc, Vec3 const& cl, Vec3 const& cu,
                                   std::function<Vec3(Vec3)> lineMaxF) {

        auto avg = [&] (Vec3 const& v1, Vec3 const& v2, Vec3 const& v3) {
            auto v_out = v1 + v2 + v3;
            return v_out / sqrt(dot(v_out, v_out));
        };

        using direc3 = std::array<Vec3, 3>;
        int ux, uy, lx, ly;
        struct Frame3 {
            direc3 d;
            direc3 p;
            Frame3(direc3 d, direc3 p ) : d(d), p(p){}
        };

        auto pcc = lineMaxF(cc),
                pcl = lineMaxF(cl),
                pcu = lineMaxF(cu);
        std::vector<Vec3> pts{ pcc, pcl, pcu };
        std::deque<Frame3> frameStack;
        frameStack.emplace_back(direc3{cc, cl, cu},
                                direc3{pcc, pcl, pcu});
        while(!frameStack.empty()) {
            Frame3 lf = frameStack.front();
            frameStack.pop_front();

// TODO fix this.
//            if (lineIntersection(lf.d[0], d1, lf.d[1], d2, lf.d[2], d3, p_ext)) {
//                if (height(lf.d[0], lf.d[1]) > eps) {
//
//                    Vec3 m_vec = avg(lf.d[0], lf.d[1], lf.d[2]);
//                    auto line_max = lineMaxF(m_vec);
//                    pts.push_back(line_max);
//                    frameStack.emplace_back(direc3{lf.d[0], m_vec, lf.d[1]},
//                                            direc3{lf.p[0], line_max, lf.p[1]});
//                    frameStack.emplace_back(direc3{lf.d[0], m_vec, lf.d[2]},
//                                            direc3{lf.p[0], line_max, lf.p[2]});
//                    frameStack.emplace_back(direc3{lf.d[1], m_vec, lf.d[2]},
//                                            direc3{lf.p[1], line_max, lf.p[2]});
//                }
//            }
        }
        return pts;
    }

    std::vector<Vec3> eps_core_set3(double eps,
                                   std::function<Vec3(Vec3)> lineMaxF) {

        auto core_set1 = eps_core_set3(eps, Vec3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{ 0.0, 1.0, 0.0 }, lineMaxF)
           , core_set2 = eps_core_set3(eps, Vec3{0.0, 0.0, 1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{ 0.0, -1.0, 0.0 }, lineMaxF)
           , core_set3 = eps_core_set3(eps, Vec3{0.0, 0.0, 1.0}, Vec3{-1.0, 0.0, 0.0}, Vec3{ 0.0, 1.0, 0.0 }, lineMaxF)
           , core_set4 = eps_core_set3(eps, Vec3{0.0, 0.0, 1.0}, Vec3{-1.0, 0.0, 0.0}, Vec3{ 0.0, -1.0, 0.0 }, lineMaxF)
           , core_set5 = eps_core_set3(eps, Vec3{0.0, 0.0, -1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{ 0.0, 1.0, 0.0 }, lineMaxF)
           , core_set6 = eps_core_set3(eps, Vec3{0.0, 0.0, -1.0}, Vec3{1.0, 0.0, 0.0}, Vec3{ 0.0, -1.0, 0.0 }, lineMaxF)
           , core_set7 = eps_core_set3(eps, Vec3{0.0, 0.0, -1.0}, Vec3{-1.0, 0.0, 0.0}, Vec3{ 0.0, 1.0, 0.0 }, lineMaxF)
           , core_set8 = eps_core_set3(eps, Vec3{0.0, 0.0, -1.0}, Vec3{-1.0, 0.0, 0.0}, Vec3{ 0.0, -1.0, 0.0 }, lineMaxF);
        core_set1.insert(core_set1.end(), core_set2.begin(), core_set2.end());
        core_set1.insert(core_set1.end(), core_set3.begin(), core_set3.end());
        core_set1.insert(core_set1.end(), core_set4.begin(), core_set4.end());
        core_set1.insert(core_set1.end(), core_set5.begin(), core_set5.end());
        core_set1.insert(core_set1.end(), core_set6.begin(), core_set6.end());
        core_set1.insert(core_set1.end(), core_set7.begin(), core_set7.end());
        core_set1.insert(core_set1.end(), core_set8.begin(), core_set8.end());
        return core_set1;
    }

  }
