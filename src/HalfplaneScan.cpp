
#include <cmath>
#include <algorithm>
#include <vector>
#include <tuple>
#include <functional>

#include "Utilities.hpp"
#include "Statistics.hpp"
#include "HalfplaneScan.hpp"


namespace pyscan {


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

    double angle(Point<> const &p1, Point<> const &p2) {
        /*
         * Finds the angle with the y-axis
         */
        double ax = get<0>(p1) - get<0>(p2);
        double ay = get<1>(p1) - get<1>(p2);
        if (ax > 0) {
            return acos(ay * invsqrt(ax * ax + ay * ay));
        } else {
            return M_PI + acos(-ay * invsqrt(ax * ax + ay * ay));
        }
    }




    template<typename F>
    std::tuple<Halfspace<>, Point<>, Point<>>
    maxHalfplane(point_it p_net_b, point_it p_net_e, point_it s_M_b, point_it s_M_e,
                 point_it s_B_b, point_it s_B_e, F func) {

        double totalM = 0;
        double totalB = 0;
        for (auto p_s = s_M_b; p_s != s_M_e; ++p_s) {
            totalM += p_s->getBlueWeight();
        }
        for (auto p_s = s_B_b; p_s != s_B_e; ++p_s) {
            totalM += p_s->getRedWeight();
        }

        if (p_net_e - p_net_b == 0) {
            return std::make_tuple(Halfspace<>(0.0, 0.0, 0.0),
                    Point<>(), Point<>());
        }
        std::vector<double> mweights(p_net_e - p_net_b - 1, 0);
        std::vector<double> bweights(p_net_e - p_net_b - 1, 0);
        std::vector<double> d_indices(p_net_e - p_net_b - 1, 0);
        std::vector<size_t> indices(p_net_e - p_net_b - 1, 0);
        Halfspace<2> maxHalfplane;
        Point<> mp1, mp2;
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
            for (auto p_s = s_M_b; p_s != s_M_e; ++p_s) {
                auto s_it = std::lower_bound(d_indices.begin(), d_indices.end(), angle(*p_b, *p_s));
                mweights[s_it - d_indices.begin()] += p_s->getBlueWeight();
            }

            for (auto p_s = s_B_b; p_s != s_B_e; ++p_s) {
                auto s_it = std::lower_bound(d_indices.begin(), d_indices.end(), angle(*p_b, *p_s));
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

    Halfspace<3> liftHalfspace(Halfspace<2> const& h2, Point<3> const& p3) {
        return {h2.fValue(), get<0>(h2), get<1>(h2),
                (1 - get<0>(p3) * get<0>(h2) - get<1>(p3) * get<1>(h2)) / get<2>(p3)};
    }

    Point<2> dropPoint(Point<3> const& fixed_point, Point<3> const& p3) {
        /*
         * Does an affine transformation from 3 to 2.
         */
        double scaling = get<2>(fixed_point) - get<2>(p3);
        return {p3.getRedWeight(), p3.getBlueWeight(),
                get<0>(p3) * get<2>(fixed_point) - get<0>(fixed_point) * get<2>(p3),
                get<1>(p3) * get<2>(fixed_point) - get<1>(fixed_point) * get<2>(p3)
        };
    }

    std::tuple<Halfspace<>, Point<>, Point<>> maxHalfplaneLin(point_it p_net_b, point_it p_net_e, point_it p_samp_b, point_it p_samp_e){
        return maxHalfplane(p_net_b, p_net_e, p_samp_b, p_samp_e, &linear);
    }

     std::tuple<Halfspace<>, Point<>, Point<>> maxHalfplaneStat(point_it p_net_b, point_it p_net_e, point_it p_samp_b, point_it p_samp_e, double rho) {
        return maxHalfplane(p_net_b, p_net_e, p_samp_b, p_samp_e, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }


    std::tuple<Halfspace<>, Point<>, Point<>> maxHalfplaneGamma(point_it p_net_b, point_it p_net_e, point_it p_samp_b, point_it p_samp_e, double rho){
        return maxHalfplane(p_net_b, p_net_e, p_samp_b, p_samp_e, [&rho](double mr, double br){
            return gamma(mr, br, rho);
        });
    }
}
