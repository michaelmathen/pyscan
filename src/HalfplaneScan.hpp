//
// Created by mmath on 9/19/17.
//

#ifndef PYSCAN_HALFPLANE_HPP
#define PYSCAN_HALFPLANE_HPP

#include <algorithm>
#include <vector>
#include <tuple>
#include <functional>
#include "Point.hpp"

namespace pyscan {


    class Halfplane {
        double a, b;
        double value;
    public:
        Halfplane(double a, double b, double val) : a(a), b(b), value(val) {}
        double getSlope() const {
            return a;
        }

        double getIntersect() const {
            return b;
        }

        double fValue() const {
            return value;
        }
    };
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

    double fastMonoAngle(Point const& p1, Point const& p2) {
        double ax = std::get<0>(p1) - std::get<0>(p2);
        double ay = std::get<1>(p1) - std::get<1>(p2);
        double f_val;
        if (ay > 0) {
            //Starts 1 is 0, 0 is 1, and -1 is 2 and then goes to .
            f_val = 1 - ax * invsqrtQuake(ax * ax + ay * ay);
        } else {
            return ax * invsqrtQuake(ax * ax + ay * ay) + 3;
        }

    }

    template <typename T, typename Compare>
    std::vector<std::size_t> sort_permutation(const std::vector<T>& vec, Compare compare) {
        std::vector<std::size_t> p(vec.size());
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(),
                  [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
        return p;
    }

    template <typename T>
    void apply_permutation_in_place(std::vector<T>& vec, std::vector<std::size_t> const& p) {
        std::vector<bool> done(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i) {
            if (done[i])
            {
                continue;
            }
            done[i] = true;
            std::size_t prev_j = i;
            std::size_t j = p[i];
            while (i != j)
            {
                std::swap(vec[prev_j], vec[j]);
                done[j] = true;
                prev_j = j;
                j = p[j];
            }
        }
    }

    template<typename F>
    Halfplane maxHalfplane(point_it p_net_b, point_it p_net_e, point_it p_samp_b, point_it p_samp_e, F func) {

        if (p_net_e - p_net_b == 0) {
            return {0, 0, 0};
        }
        std::vector<double> weights(p_net_e - p_net_b - 1, 0);
        std::vector<double> d_indices(p_net_e - p_net_b - 1, 0);
        std::vector<size_t> indices(p_net_e - p_net_b - 1, 0);
        for (auto p_b = p_net_b; p_b != p_net_e; p_b++) {

            auto dind_it = d_indices.begin();
            auto ind_it =  indices.begin();
            auto p_b_i = p_net_b;
            while (p_b_i != p_net_e) {
                if (p_b_i != p_b) {
                    *ind_it = p_b_i - p_net_b;
                    *dind_it = fastMonoAngle(*p_b, *p_b_i);
                    ++ind_it;
                    ++dind_it;
                }
                ++p_b_i;
            }
            auto perm = sort_permutation(d_indices.begin(), d_indices.end(), std::less);
            // Now need to sort indices and weights.
            //Apply to indices and weights

        }
    }
}
#endif //PYSCAN_HALFPLANE_HPP
