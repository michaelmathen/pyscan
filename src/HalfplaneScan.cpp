
#include <cmath>
#include <algorithm>
#include <vector>
#include <tuple>
#include <functional>

#include "Statistics.hpp"
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

    double angle(Point<double> const &p1, Point<double> const &p2) {
        /*
         * Finds the angle with the y-axis
         */
        double ax = get<0>(p1) - get<0>(p2);
        double ay = get<1>(p1) - get<1>(p2);
        if (ax > 0) {
            return acos(ay * invsqrtQuake(ax * ax + ay * ay));
        } else {
            return M_PI + acos(-ay * invsqrtQuake(ax * ax + ay * ay));
        }
    }




    template<typename F>
    std::tuple<Halfspace<>, Point<double>, Point<double>>
    maxHalfplane(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e, F func) {

        double totalM = 0;
        double totalB = 0;
        for (auto p_s = p_samp_b; p_s != p_samp_e; ++p_s) {
            totalM += p_s->getBlueWeight();
            totalB += p_s->getRedWeight();
        }

        if (p_net_e - p_net_b == 0) {
            return std::make_tuple(Halfspace<>(0.0, 0.0, 0.0),
                    Point<double>(), Point<double>());
        }
        std::vector<double> mweights(p_net_e - p_net_b - 1, 0);
        std::vector<double> bweights(p_net_e - p_net_b - 1, 0);
        std::vector<double> d_indices(p_net_e - p_net_b - 1, 0);
        std::vector<size_t> indices(p_net_e - p_net_b - 1, 0);
        Halfspace<2> maxHalfplane;
        Point<double> mp1, mp2;
        for (auto p_b = p_net_b; p_b != p_net_e; p_b++) {

            auto dind_it = d_indices.begin();
            auto ind_it =  indices.begin();
            size_t it = 0;
            auto p_b_i = p_net_b;
            while (p_b_i != p_net_e) {
                if (p_b_i != p_b) {
                    *ind_it = it;
                    *dind_it = angle(*p_b, *p_b_i);
                    ++ind_it;
                    ++dind_it;
                }
                ++p_b_i;
                ++it;
            }
            //Sort by the angles
            auto perm = sort_permutation(d_indices, std::less<double>());

            // Now need to sort indices and weights.
            //Apply to indices and weights
            apply_permutation_in_place(indices, perm);

            //Now update all the weights as the angle changes over all the points
            for (auto p_s = p_samp_b; p_s != p_samp_e; ++p_s) {
                auto s_it = std::lower_bound(d_indices.begin(), d_indices.end(), angle(*p_b, *p_s));
                mweights[s_it - d_indices.begin()] += p_s->getBlueWeight();
                bweights[s_it - d_indices.begin()] += p_s->getRedWeight();
            }

            // Now find the value of the first halfspace.
            double m_cum = 0, b_cum = 0;
            auto halfplane_end_it = std::lower_bound(d_indices.begin(), d_indices.end(), M_PI);
            {
                auto mW_e = mweights.begin() + (halfplane_end_it - d_indices.begin());
                auto bW_e = bweights.begin() + (halfplane_end_it - d_indices.begin());
                for (auto h_e_it = halfplane_end_it; h_e_it != d_indices.end();) {
                    m_cum += *mW_e;
                    b_cum += *bW_e;
                    ++mW_e, ++bW_e, ++h_e_it;
                }
            }
            // Now do the scan.
            {
                double maxF = func(m_cum / totalM, b_cum / totalB);
                decltype(indices.begin()) curr_it;
                decltype(indices.begin()) max_it;
                if (*halfplane_end_it < *d_indices.begin() + M_PI) {
                    curr_it = indices.begin() + (halfplane_end_it - d_indices.begin());
                } else {
                    curr_it = indices.begin();
                }
                max_it = curr_it;
                auto b_d_ind = d_indices.begin();
                auto e_d_ind = halfplane_end_it;
                auto mW_b = mweights.begin();
                auto bW_b = bweights.begin();
                auto mW_e = mweights.begin() + (halfplane_end_it - d_indices.begin());
                auto bW_e = bweights.begin() + (halfplane_end_it - d_indices.begin());

                while (b_d_ind != halfplane_end_it || e_d_ind != d_indices.end()) {

                    if (b_d_ind == halfplane_end_it || *b_d_ind + M_PI > *e_d_ind) {
                        m_cum -= *mW_e;
                        b_cum -= *bW_e;
                        curr_it = indices.begin() + (e_d_ind - d_indices.begin());
                        ++e_d_ind, ++mW_e, ++bW_e;

                    } else if (e_d_ind == d_indices.end() || *b_d_ind + M_PI < *e_d_ind){
                        m_cum += *mW_b;
                        b_cum += *bW_b;
                        curr_it = indices.begin() + (b_d_ind - d_indices.begin());
                        ++b_d_ind, ++mW_b, ++bW_b;
                    } else {
                        m_cum += *mW_b;
                        b_cum += *bW_b;
                        m_cum -= *mW_e;
                        b_cum -= *bW_e;
                        curr_it = indices.begin() + (b_d_ind - d_indices.begin());
                        ++b_d_ind, ++mW_b, ++bW_b,
                                ++e_d_ind, ++mW_e, ++bW_e;
                    }
                    double tmpF = func(m_cum / totalM, b_cum / totalB);
                    if (tmpF >= maxF) {
                        maxF = tmpF;
                        max_it = curr_it;
                    }
                }
                if (maxF >= maxHalfplane.fValue()) {
                    mp1 =  *p_b;
                    mp2 = *(p_b + *max_it);
                    double a = (get<1>(mp2) - get<1>(mp1)) / (get<0>(mp1) * get<1>(mp2) - get<1>(mp1) * get<0>(mp2));
                    double b = (1 - a / get<0>(mp1)) / get<0>(mp2);
                    maxHalfplane = Halfspace<>(a, b, maxF);
                }
                //Now figure out which point the max angle corresponds to.
            }
        }
        return std::make_tuple(maxHalfplane, mp1, mp2);
    }


    std::tuple<Halfspace<>, Point<double>, Point<double>> maxHalfplaneLin(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e){
        return maxHalfplane(p_net_b, p_net_e, p_samp_b, p_samp_e, &linear);
    }

     std::tuple<Halfspace<>, Point<double>, Point<double>> maxHalfplaneStat(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e, double rho) {
        return maxHalfplane(p_net_b, p_net_e, p_samp_b, p_samp_e, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }


    std::tuple<Halfspace<>, Point<double>, Point<double>> maxHalfplaneGamma(point_d_it p_net_b, point_d_it p_net_e, point_d_it p_samp_b, point_d_it p_samp_e, double rho){
        return maxHalfplane(p_net_b, p_net_e, p_samp_b, p_samp_e, [&rho](double mr, double br){
            return gamma(mr, br, rho);
        });
    }
}