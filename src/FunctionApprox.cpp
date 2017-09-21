//
// Created by mmath on 7/7/17.
//

#include <cmath>
#include <deque>
#include <algorithm>
#include <functional>
#include <tuple>


inline void dot(double A[4], double B[2], double X[2]){
    X[0] = A[0] * B[0] + A[1] * B[1];
    X[1] = A[2] * B[0] + A[3] * B[1];
}


inline bool lineIntersection(double alpha, double theta, double mi, double bi,
                             double mj, double bj, double& m, double& b) {
    /*
    Find the point where lines
    l(t) = (sin(alpha), -cos(alpha))t + (mi, bi)
    l(t) = (sin(theta + alpha), -cost(theta + alpha))t  + (mj, bj)
    intersect. This point corresponds to the maximum possible point
    in a convex set.
     */
    if (theta <= 0) {
        return true;
    }

    //	 double A1[4] = {-sin(alpha), -sin(theta + alpha),
    //			 	 	 	 	 	 	 cos(alpha),  cos(theta + alpha) };
    //   double A2[4] = {cos(theta + alpha), sin(theta + alpha),
    //	                -cos(alpha)       , -sin(alpha)};
    double B[2] = {bi - bj, mi - mj};
    double A[4] = {-sin(2*alpha + theta), -2 * sin(alpha) * sin(alpha + theta),
                   2 * cos(alpha) * cos(theta + alpha), sin(2 * alpha + theta)
    };
    double X[2];
    dot(A, B, X);

    double det = sin(theta);
    //double t_j = (cos(theta + alpha) * (bi - bj) + sin(theta + alpha) *(mi - mj)) / sin(theta);
    //double mdiff = mi - mj;
    //double bdiff = bi - bj;
    //cout << -sin(alpha) * mdiff - cos(alpha) * bdiff << endl;
    //double t_i = (-sin(alpha) * mdiff - cos(alpha) * bdiff) / sin(theta);
    //cout << t_i << endl;
    //cout << sin(alpha + theta) << endl;
    //b =  t_i * sin(alpha + theta) + bi;
    //m =  -t_i * cos(alpha + theta) + mi;
    b = 1 / 2.0 * (1 / det * X[0] + bi + bj);
    m = 1 / 2.0 * (1 / det * X[1] + mi + mj);
    return false;
}




double approximateHull(double mi, double bi,
                       double mj, double bj,
                       double alpha, double theta, double eps,
                       std::function<double(double, double)> phi, //function to maximize
                       std::function<std::tuple<double, double>(double)> lineMaxF)
{
    //int ux, uy, lx, ly;
    struct Frame {
        double alpha, theta, mi, bi, mj, bj;
        Frame(double a, double t, double m_i, double b_i, double m_j, double b_j) :
                alpha(a), theta(t), mi(m_i), bi(b_i), mj(m_j), bj(b_j) {}
    };
    double maxRValue = 0;

    std::deque<Frame> frameStack;
    // TODO double check debug to see if there is an issue with an infinite singularity here. Might need to change start.
    frameStack.push_back(Frame(alpha, theta, mi, bi, mj, bj));

    while(!frameStack.empty()) {
        Frame lf = frameStack.front();
        frameStack.pop_front();
        //cout << lf.alpha << " " << lf.theta << " " << lf.bi << " " << lf.mi << " " << lf.bj << " " << lf.mj << endl;
        double mExt, bExt;
        if (lineIntersection(lf.alpha, lf.theta, lf.mi, lf.bi, lf.mj, lf.bj, mExt, bExt)) {
            continue;
        }
        double vw = phi(mExt, bExt);
        double vi = phi(lf.mi, lf.bi);
        double vj = phi(lf.mj, lf.bj);
// 		 cout << "i phi(" << lf.bi << " " << lf.mi << ")=" << vi << endl;
// 		 cout << "u phi(" << bExt << " " << mExt << ")=" << vw << endl;
// 		 cout << "j phi(" << lf.bj << " " << lf.mj << ")=" << vj << endl;
        if (std::max({vi, vj, maxRValue}) + eps  > vw) {
            //cout << "discarding result" << endl;
            // No value in this triangle is large enough to be interesting
            continue;
        } else {
            //This triangle is worth evaluating
            auto lineMax = lineMaxF(lf.alpha + lf.theta / 2);
            frameStack.push_back(Frame(lf.alpha , lf.theta/2,
                                       std::get<0>(lineMax), std::get<1>(lineMax), lf.mj, lf.bj));
            frameStack.push_back(Frame(lf.alpha + lf.theta / 2, lf.theta / 2 ,
                                       lf.mi, lf.bi, std::get<0>(lineMax), std::get<1>(lineMax)));
            //double m = get<0>(lineMax);
            //double b = get<1>(lineMax);
            //cout << "n phi(" << b << " " << m << ")=" << phi(m, b) << endl;
        }
    }
    return maxRValue;
}
