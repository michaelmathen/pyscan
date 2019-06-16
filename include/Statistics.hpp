//
// Created by mmath on 7/7/17.
//

#ifndef PYSCAN_STATISTICS_HPP
#define PYSCAN_STATISTICS_HPP
#include <cmath>
#include <limits>
#include <iostream>


namespace pyscan {

    inline double kulldorff(double mr, double br, double rho) {
        if (mr < rho || br < rho || br > 1 - rho || mr > 1 - rho) {
            return 0;
        }
        if (fabs(1 - fabs(mr / br)) <= std::numeric_limits<double>::epsilon()) {
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

    inline double regularized_kulldorff(double mr, double br, double rho) {
        if (mr == 0) {
            return log(1 / (1 - br + rho));
        } else if (mr == 1) {
            return log(1 / (br + rho));
        } else {
            return mr * log(mr / (br + rho)) + (1 - mr) * log((1 - mr) / (1 - br + rho));
        }
    }


    inline double gamma(double mr, double br, double rho) {
        if (br < rho || br > 1 - rho)
            return 0;
        if (br <= 0 || br >= 1)
            return std::numeric_limits<double>::infinity();
        return (mr - br) * (mr - br) / (br * (1 - br));
    }

    inline double bound(double value, double rho) {
        return std::min(1 - rho, std::max(rho, value));
    }


    inline discrepancy_func_t get_bernoulli(double rho) {
        //G = 1.0;
        auto bernoulli = [rho] (double m_sub, double m_total, double b_sub, double b_total) {

            double p = m_sub / (b_sub + m_sub);
            double q = (m_total - m_sub) / (b_total - b_sub + m_total - m_sub);
            return m_sub * log(p)
                   + b_sub * log(1 - p)
                   + (m_total - m_sub) * log(q)
                   + (b_total - b_sub) * log(1 - q)
                   - m_total * log(m_total / (b_total + m_total))
                   - b_total * log(1 - m_total/ (m_total + b_total));
        };
        discrepancy_func_t b_f = bernoulli;
        return b_f;
    }

    inline discrepancy_func_t get_bernoulli_single_sample(double rho) {
        //G = 1.0;
        auto bernoulli = [rho] (double m_sub, double m_total, double b_sub, double b_total) {
            //Bound the function.
            // This changes the behaviour when the region is very large or small to prevent the function's lipshitz
            // parameter from blowing up.
//            return m_sub * log(m_sub / m_total) + (m_total - m_sub) * log(1 - m_sub/m_total) +
//                    b_sub * log(b_sub / m_sub) + (b_total - b_sub) * log(1 - b_sub / b_total)
//                    - m_total * log(m_sub / m_total) - b_total * log(b_sub / b_total);
            //mu(Z) = b_sub
            //n_z = m_sub
            //n_G = m_total
            //mu(G) = b_total
            //double m_r = m_sub / m_total;
            //double b_r = b_sub / b_total;
            double p = m_sub / (b_sub);
            double q = (m_total - m_sub) / (b_total - b_sub);

            return m_sub * log(p)
                   + (b_sub - m_sub) * log(1 - p)
                   + (m_total - m_sub) * log(q)
                   + (b_total - m_total -(b_sub - m_sub)) * log(1 - q)
                   - m_total * log(m_total / b_total)
                   - (b_total - m_total) * log(1 - m_total / b_total);
        };
        discrepancy_func_t b_f = bernoulli;
        return b_f;
    }

    inline double linear(double mr, double br) {
        return  fabs(mr - br);
    }

    // template<typename Reg, typename F>
    // double evaluateRegion(std::vector<LPoint<>> const& m_pts, std::vector<LPoint<>> const& b_pts, Reg const& reg, F func) {
    //   double m_total = computeLabelTotal(m_pts.begin(), m_pts.end(), getMeasured);
    //   double b_total = computeLabelTotal(b_pts.begin(), b_pts.end(), getBaseline);
    //   auto filterF = [&] (Point<> const& pt) {
    //     return reg.contains(pt);
    //   };
    //   double m_curr = computeLabelTotalF(m_pts.begin(), m_pts.end(), getMeasured,
    //                                     filterF);
    //   double b_curr = computeLabelTotalF(b_pts.begin(), b_pts.end(), getBaseline,
    //                                     filterF);
    //   return func(m_curr / m_total, b_curr / b_total);
    // }

    // template<typename Reg, typename F>
    // double evaluateRegion(std::vector<Point<>>& m_pts, std::vector<Point<>>& b_pts, Reg const& reg, F func) {
    //     double m_curr = 0;
    //     double m_total = 0;
    //     for (auto p = m_pts.begin(); p != m_pts.end(); p++) {
    //         if (reg.contains(*p)) {
    //             m_curr += getMeasured(*p);
    //         }
    //         m_total += getMeasured(*p);
    //     }

    //     double b_curr = 0;
    //     double b_total = 0;
    //     for (auto p = b_pts.begin(); p != b_pts.end(); p++) {
    //         if (reg.contains(*p)) {
    //             b_curr += getBaseline(*p);
    //         }
    //         b_total += getBaseline(*p);
    //     }
    //     return func(m_curr / m_total, b_curr / b_total);
    // }
}
#endif //PYSCAN_STATISTICS_HPP
