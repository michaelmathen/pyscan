
#include <cmath>
#include "HalfplaneScan.hpp"

namespace pyscan {


    double invsqrtQuake( double number ) {
        double y = number;
        double x2 = y * 0.5;
        std::int64_t i = *(std::int64_t *) &y;
        // The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
        i = 0x5fe6eb50c7b537a9 - (i >> 1);
        y = *(double *) &i;
        y = y * (1.5 - (x2 * y * y));   // 1st iteration
        y  = y * ( 1.5 - ( x2 * y * y ) );   // 2nd iteration, this can be removed
        return y;
    }

    double fastMonoAngle(Point<double> const& p1, Point<double> const& p2) {
        double ax = p1.getX() - p2.getX();
        double ay = p1.getY() - p2.getY();
        double f_val;
        if (ay > 0) {
            //Starts 1 is 0, 0 is 1, and -1 is 2 and then goes to .
            f_val = 1 - ax * invsqrtQuake(ax * ax + ay * ay);
        } else {
            return ax * invsqrtQuake(ax * ax + ay * ay) + 3;
        }

    }

    double angle(Point<double> const &p1, Point<double> const &p2) {
        /*
         * Finds the angle with the y-axis
         */
        double ax = p1.getX() - p2.getX();
        double ay = p1.getY() - p2.getY();
        if (ax > 0) {
            return acos(ay * invsqrtQuake(ax * ax + ay * ay));
        } else {
            return M_PI + acos(-ay * invsqrtQuake(ax * ax + ay * ay));
        }
    }

}