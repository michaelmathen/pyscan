#include "Point.hpp"

namespace pyscan {


    double getMeasured(Point<> const& pt) {
        return pt.getRedWeight();
    }

    double getBaseline(Point<> const& pt) {
        return pt.getBlueWeight();
    }


    void getLoc(Point<> const& pt, double& x1, double& x2) {
        x1 = get<0>(pt);
        x2 = get<1>(pt);
    }



    double getX(Point<> const& pt) {
        return get<0>(pt);
    }

    double getY(Point<> const& pt) {
        return get<1>(pt);
    }

    bool colinear(Point<> const& pt1,
                  Point<> const& pt2,
                  Point<> const& pt3){
        double x1, x2, x3, y1, y2, y3;
        getLoc(pt1, x1, y1);
        getLoc(pt2, x2, y2);
        getLoc(pt3, x3, y3);

        double
                a11 = x2 - x1, a12 = y2 - y1,
                a21 = x2 - x3, a22 = y2 - y3;
        return (a11 * a22 - a12 * a21 == 0);
    }
}
