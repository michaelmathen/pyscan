//
// Created by mmath on 7/7/17.
//

#ifndef PYSCAN_STATISTICS_HPP
#define PYSCAN_STATISTICS_HPP
#include <cmath>
#include <limits>
#include <iostream>

#include "DiskScan.hpp"
#include "Point.hpp"

namespace pyscan {

    inline double kulldorff(double mr, double br, double rho) {
        if (mr < rho || br < rho || br > 1 - rho || mr > 1 - rho) {
            return 0;
        }
        if (abs(1 - abs(mr / br)) <= std::numeric_limits<double>::epsilon()) {
            return 0;
        }
        if (mr == 0) {
            if (br == 1) {
                return std::numeric_limits<double>::infinity();
            } else {
                return log(1 / (1 - br));
            }
        }
        if (mr == 1) {
            if (br == 0) {
                return std::numeric_limits<double>::infinity();
            } else {
                return log(1 / br);
            }
        }
        if (br == 0 || br == 1) {
            return std::numeric_limits<double>::infinity();
        }

        return mr * log(mr / br) + (1 - mr) * log((1 - mr) / (1 - br));
    }


    inline double gamma(double mr, double br, double rho) {
        if (br < rho || br > 1 - rho)
            return 0;
        if (br <= 0 || br >= 1)
            return std::numeric_limits<double>::infinity();
        return (mr - br) * (mr - br) / (br * (1 - br));
    }

    inline double linear(double mr, double br) {
        return  abs(mr - br);
    }

    template<typename Reg, typename F>
    double evaluateRegion(std::vector<LPoint<>> const& m_pts, std::vector<LPoint<>> const& b_pts, Reg const& reg, F func) {
      double m_total = computeLabelTotal(m_pts.begin(), m_pts.end(), getMeasured);
      double b_total = computeLabelTotal(b_pts.begin(), b_pts.end(), getBaseline);
      auto filterF = [&] (Point<> const& pt) {
        return reg.contains(pt);
      };
      double m_curr = computeLabelTotalF(m_pts.begin(), m_pts.end(), getMeasured,
                                        filterF);
      double b_curr = computeLabelTotalF(b_pts.begin(), b_pts.end(), getBaseline,
                                        filterF);
      return func(m_curr / m_total, b_curr / b_total);
    }

    template<typename Reg, typename F>
    double evaluateRegion(std::vector<Point<>>& m_pts, std::vector<Point<>>& b_pts, Reg const& reg, F func) {
        double m_curr = 0;
        double m_total = 0;
        for (auto p = m_pts.begin(); p != m_pts.end(); p++) {
            if (reg.contains(*p)) {
                m_curr += getMeasured(*p);
            }
            m_total += getMeasured(*p);
        }

        double b_curr = 0;
        double b_total = 0;
        for (auto p = b_pts.begin(); p != b_pts.end(); p++) {
            if (reg.contains(*p)) {
                b_curr += getBaseline(*p);
            }
            b_total += getBaseline(*p);
        }
        return func(m_curr / m_total, b_curr / b_total);
    }
}
#endif //PYSCAN_STATISTICS_HPP
