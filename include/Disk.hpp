#ifndef __DISK_H__
#define __DISK_H__

#include "Range.hpp"
#include "Point.hpp"

namespace pyscan {

    class Disk : public Range<2> {
    public:
        Disk(double x_, double y_, double r_)
                : origin(x_, y_, 1.0), radius(r_) {}

        Disk()
                : origin(0.0, 0.0, 1.0), radius(0.0) {}

        Disk(const pt2_t &pt1, const pt2_t &pt2, const pt2_t &pt3) {
            double x1, x2, x3, y1, y2, y3;
            x1 = pt1(0);
            y1 = pt1(1);
            x2 = pt2(0);
            y2 = pt2(1);
            x3 = pt3(0);
            y3 = pt3(1);

            double a11 = x2 - x1;
            double a12 = y2 - y1;
            double a21 = x2 - x3;
            double a22 = y2 - y3;
            double b1 = (y2 * y2 + x2 * x2 - y1 * y1 - x1 * x1) / 2.0;
            double b2 = (y2 * y2 + x2 * x2 - y3 * y3 - x3 * x3) / 2.0;

            double detA = util::det2(a11, a12, a21, a22);
            // inverse of A
            double ai11 = a22 / detA;
            double ai12 = -a12 / detA;
            double ai21 = -a21 / detA;
            double ai22 = a11 / detA;

            // A^{-1} b = x
            double a, b;
            a = ai11 * b1 + ai12 * b2;
            b = ai21 * b1 + ai22 * b2;
            origin = pt2_t(a, b, 1.0);
            radius = sqrt((x1 - a) * (x1 - a) + (y1 - b) * (y1 - b));
        }

        Disk(const pt2_t &p1, const pt2_t &p2, double r_) : origin(0.0, 0.0, 1.0), radius(r_) {
            double h = 1 / 2.0 * sqrt(4 * r_ * r_ - p1.square_dist(p2));
            pt2_t midpoint = 1 / 2.0 * (p1 + p2);
            pt2_t par = p1.direction(p2);
            double px = -par(1);
            double py = par(0);
            origin = Point<2>(midpoint(0) + h * px, midpoint(1) + h * py, 1.0);
        }

        inline double getRadius() const {
            return radius;
        }

        inline pt2_t getOrigin() const {
            return origin;
        }

        virtual std::string str() const {
            std::stringstream ss;
            ss << "Disk(";
            ss << this->getOrigin();
            ss << ", ";
            ss << this->getRadius();
            ss << ")";
            return ss.str();
        }

        inline bool contains(const pt2_t &pt) const final {
            return util::alte(pt.square_dist(origin), radius * radius);
        }

        inline bool intersects_segment(const pt2_t &p1, const pt2_t &p2) const final {
            //Checks to see if this intersects the segment.
            return util::alte(origin.square_dist(p1, p2) , radius * radius);
        }

    private:
        pt2_t origin;
        double radius;
    };

}

#endif
