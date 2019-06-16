/*
 * Created by Michael Matheny on 4/25/19.
 * at the University of Utah
 * email: michaelmathen@gmail.com
 * website: https://mmath.dev/
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <tuple>
#include <functional>
#include <memory>
#include <gsl/gsl_multimin.h>

#include "Sampling.hpp"
#include "KernelScanning.hpp"
#include "Gridding.hpp"
#include "Utilities.hpp"


namespace pyscan {



//    kernel_func_t rectangle_kernel(double deviation) {
//        return [deviation] (double dist) {
//            return ;
//        };
//    }
//
//    kernel_func_t triangle_kernel(double deviation);
//
//    kernel_func_t epanechnikov_kernel(double deviation);
//


    double
    unwrap_f(const gsl_vector *v, void *params)
    {
        auto fp = static_cast<Bernoulli_Disk*>(params);

        double p = gsl_vector_get(v, 0);
        double q = gsl_vector_get(v, 1);
        return -fp->alternative_hyp(p, q);
    }

    /* The gradient of f, df = (df/dx, df/dy). */
    void
    unwrap_df (const gsl_vector *v, void *params,
           gsl_vector *df)
    {
        auto fp = static_cast<Bernoulli_Disk*>(params);

        double p = gsl_vector_get(v, 0);
        double q = gsl_vector_get(v, 1);

        auto [dfp, dfq] = fp->diff(p, q);
        gsl_vector_set(df, 0, -dfp);
        gsl_vector_set(df, 1, -dfq);
    }

    /* Compute both f and df together. */
    void
    unwrap_fdf (const gsl_vector *x, void *params,
            double *f, gsl_vector *df)
    {
        *f = unwrap_f(x, params);
        unwrap_df(x, params, df);
    }

    std::tuple<double, double, double> find_pq_poi(
            double p_init,
            double q_init,
            const Bernoulli_Disk& disc_f) {
        size_t iter = 0;
        int status = GSL_CONTINUE;

        const gsl_multimin_fdfminimizer_type *T;
        gsl_multimin_fdfminimizer *s;

        gsl_vector *x;
        gsl_multimin_function_fdf my_func;

        my_func.n = 2;
        my_func.f = unwrap_f;
        my_func.df = unwrap_df;
        my_func.fdf = unwrap_fdf;
        my_func.params = (void*)&disc_f;

        /*Initialize the vector x to be p_init and q_init*/
        x = gsl_vector_alloc (2);
        gsl_vector_set (x, 0, p_init);
        gsl_vector_set (x, 1, q_init);

        T = gsl_multimin_fdfminimizer_conjugate_fr;
        s = gsl_multimin_fdfminimizer_alloc (T, 2);

        /*
         * The initial step-size is chosen as 0.01, and the line minimization parameter is set at 0.0001.
         */
        gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

        double p_o = p_init, q_o = q_init;
        do
        {
            iter++;

            status = gsl_multimin_fdfminimizer_iterate (s);

            if (status)
                break;
            /*
             * The program terminates when the norm of the gradient has been reduced below 0.001.
             */
            double p = gsl_vector_get(s->x, 0);
            double q = gsl_vector_get(s->x, 1);
            bool trig = true;
            if (p <= 0) {
                p = 1e-4;
            } else if (p > 1) {
                p = 1 - 1e-4;
            } else if (q <= 0) {
                q = 1e-4;
            } else if (q >= 1) {
                q = 1 - 1e-4;
            } else if (std::isnan(q) || std::isnan(p)) {
                //restart at the last location.
                p = p_o;
                q = q_o;
            } else {
                trig = false;
            }

            if (trig) {
                gsl_vector_set (x, 0, p);
                gsl_vector_set (x, 1, q);
                gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);
            } else {
                status = gsl_multimin_test_gradient(s->gradient, 1e-4);
            }
            p_o = p;
            q_o = q;
        } while (status == GSL_CONTINUE && iter < 100);

        double p = gsl_vector_get(s->x, 0);
        double q = gsl_vector_get(s->x, 1);
        gsl_multimin_fdfminimizer_free (s);
        gsl_vector_free (x);
        double f_val = disc_f.lrt(p, q);
        //If this region does not have a higher than normal measured than baseline points
        return std::make_tuple(p, q, f_val);
    }

    std::tuple<double, double, double> measure_kernel(
            const pt2_t& center,
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double bandwidth) {

        double red_tot = computeTotal(measured);
        double blue_tot = computeTotal(baseline);
        auto kern = [](double dist, double bandwidth) {
            return exp(-dist * dist / (bandwidth * bandwidth));
        };
        Bernoulli_Disk disc(red_tot, blue_tot, bandwidth, kern);

        std::vector<double> m_weights;
        std::vector<double> b_weights;
        std::vector<double> m_radii;
        std::vector<double> b_radii;

        for (auto& p : measured) {
            m_weights.emplace_back(p.get_weight());
            m_radii.emplace_back(center.dist(p));
        }
        for (auto& p : baseline) {
            b_weights.emplace_back(p.get_weight());
            b_radii.emplace_back(center.dist(p));
        }
        disc.set_weights(m_weights, b_weights);
        disc.set_radii(m_radii, b_radii);

        return find_pq_poi(.6, .5, disc);
    }

    struct Kernel {
        double operator()(double dist, double bandwidth) {
            return exp(-dist * dist / (bandwidth * bandwidth));
        }
    };


    std::tuple<Disk, double> max_kernel_slow_internal(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double disk_r,
            Bernoulli_Disk& disc,
            bbox_t const& full_bb) {


        Disk curr_max;
        double max_stat = 0.0;

        std::vector<double> m_weights;
        std::vector<double> b_weights;
        for (auto& p : measured) {
            m_weights.emplace_back(p.get_weight());
        }
        for (auto& p : baseline) {
            b_weights.emplace_back(p.get_weight());
        }
        disc.set_weights(m_weights, b_weights);
        std::vector<double> m_radii(measured.size(), 0.0);
        std::vector<double> b_radii(baseline.size(), 0.0);

        double p_init = .6;
        double q_init = .5;

        auto [mnx, mny, mxx, mxy] = full_bb;
        for (double x = mnx; x < mxx; ) {
            for (double y = mny; y < mxy; ) {

                pt2_t center(x, y, 1.0);
                Disk disk(x, y, disk_r);
                for (size_t i = 0; i < measured.size(); i++) {
                    m_radii[i] = center.dist(measured[i]);
                }
                for (size_t i = 0; i < baseline.size(); i++) {
                    b_radii[i] = center.dist(baseline[i]);
                }

                disc.set_radii(m_radii, b_radii);

                auto[p, q, fval] = find_pq_poi(p_init, q_init, disc);


                if (std::isnan(p) || std::isnan(q) || p <= 0 || q <= 0 || p >= 1 || q >= 1) {
                    p_init = .6;
                    q_init = .5;
                } else {
                    p_init = p;
                    q_init = q;
                }
                if (max_stat < fval) {
                    max_stat = fval;
                    curr_max = disk;
                }
                y = y + grid_res;
            }
            x = x + grid_res;
        }
        return std::make_tuple(curr_max, max_stat);
    }

    std::tuple<Disk, double> max_kernel_slow(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double disk_r,
            double bandwidth) {
        double red_tot = computeTotal(measured);
        double blue_tot = computeTotal(baseline);
        Kernel kern;
        Bernoulli_Disk disc(red_tot, blue_tot, bandwidth, kern);
        auto bb_op = bbox(measured, baseline);
        Disk curr_max;
        double max_stat = 0.0;

        if (!bb_op.has_value()) {
            return std::make_tuple(curr_max, max_stat);
        }
        auto full_bb = bb_op.value();
        auto [mnx, mny, mxx, mxy] = full_bb;
        auto edited_bb = std::make_tuple(mnx - disk_r, mny - disk_r, mxx + disk_r, mxy + disk_r);
        return max_kernel_slow_internal(measured, baseline, grid_res, disk_r, disc, edited_bb);
    }


    template<bool enable_prune, bool adaptive_grid>
    std::tuple<Disk, double> max_kernel_internal(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth) {


        double red_tot = computeTotal(measured);
        double blue_tot = computeTotal(baseline);
        Kernel kern;
        Bernoulli_Disk disc(red_tot, blue_tot, bandwidth, kern);

        Disk curr_max;
        double max_stat = 0.0;

        auto bb_op = bbox(measured, baseline);
        if (!bb_op.has_value()) {
            return std::make_tuple(curr_max, max_stat);
        }
        auto full_bb = bb_op.value();

        for (size_t x_off = 0; x_off < 3; x_off++) {
            for (size_t y_off = 0; y_off < 3; y_off++) {
                auto bb = std::make_tuple(
                        std::get<0>(full_bb) + x_off * radius_size - 2 * radius_size,
                        std::get<1>(full_bb) + y_off * radius_size - 2 * radius_size,
                        std::get<2>(full_bb) + x_off * radius_size,
                        std::get<3>(full_bb) + y_off * radius_size);
                //Define a coarse grid.
                SparseGrid<wpt2_t> grid_red(bb, measured, radius_size * 3), grid_blue(bb, baseline, radius_size * 3);
                std::vector<uint64_t> keys;
                {
                    for (auto b = grid_red.begin(); b != grid_red.end(); b++) {
                        keys.emplace_back(b->first);
                    }
                    for (auto b = grid_blue.begin(); b != grid_blue.end(); b++) {
                        keys.emplace_back(b->first);
                    }
                    std::sort(keys.begin(), keys.end());
                    auto nend = std::unique(keys.begin(), keys.end());
                    keys.erase(nend, keys.end());
                }

                //Consider cells with points in them, so they must have $\eps n$ values
                for (auto k : keys) {

                    auto[lx, ly] = grid_red.get_lower_corner(k);
                    auto new_bb = std::make_tuple(
                            lx + radius_size,
                            ly + radius_size,
                            lx + 2 * radius_size,
                            ly + 2 * radius_size);
                    //std::cout << lx + radius_size << " " << ly + radius_size << " " << lx + 2* radius_size << " " << ly + 2 * radius_size << std::endl;

                    wpoint_list_t measured_local;
                    wpoint_list_t baseline_local;
                    auto[r_it, r_e_it] = grid_red(k);

                    double subgrid_red_total = 0.0;
                    for (auto b = r_it; b != r_e_it; b++) {
                        subgrid_red_total += b->second.get_weight();
                        measured_local.emplace_back(b->second);
                    }
                    auto[b_it, b_e_it] = grid_blue(k);
                    double subgrid_blue_total = 0.0;
                    for (auto b = b_it; b != b_e_it; b++) {
                        subgrid_blue_total += b->second.get_weight();
                        baseline_local.emplace_back(b->second);
                    }
                    double dense_grid_res;
                    if (adaptive_grid) {
                        dense_grid_res = grid_res * (red_tot + blue_tot) / (subgrid_blue_total + subgrid_red_total);
                    } else {
                        dense_grid_res = grid_res;
                    }
                    Disk max_d;
                    double max_v;
                    if (enable_prune) {
                        std::tie(max_d, max_v) = max_kernel_slow_internal(
                                measured_local,
                                baseline_local,
                                dense_grid_res,
                                radius_size,
                                disc, new_bb);
                    } else {
                        std::tie(max_d, max_v) = max_kernel_slow_internal(
                                measured,
                                baseline,
                                dense_grid_res,
                                radius_size,
                                disc, new_bb);

                    }
                    if (max_v > max_stat) {
                        max_stat = max_v;
                        curr_max = max_d;
                    }
                }
            }
        }

        return std::make_tuple(curr_max, max_stat);
    }

    std::tuple<Disk, double> max_kernel_slow2(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth) {
        return max_kernel_internal<false, false>(measured, baseline, grid_res, radius_size, bandwidth);
    }

     std::tuple<Disk, double> max_kernel_adaptive(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth) {
        return max_kernel_internal<false, true>(measured, baseline, grid_res, radius_size, bandwidth);
     }

     std::tuple<Disk, double> max_kernel_prune_far(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth) {
        return max_kernel_internal<true, false>(measured, baseline, grid_res, radius_size, bandwidth);
    }

    std::tuple<Disk, double> max_kernel(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth) {
        return max_kernel_internal<true, true>(measured, baseline, grid_res, radius_size, bandwidth);
    }

    void max_kernel_slow_centers(
            point_list_t & centers,
            double grid_res,
            bbox_t const& full_bb) {

        auto [mnx, mny, mxx, mxy] = full_bb;
        for (double x = mnx; x < mxx; ) {
            for (double y = mny; y < mxy; ) {
                pt2_t center(x, y, 1.0);
                centers.emplace_back(center);
                y = y + grid_res;
            }
            x = x + grid_res;
        }
    }

    point_list_t kernel_centers_approximate(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth) {


        double red_tot = computeTotal(measured);
        double blue_tot = computeTotal(baseline);

        auto bb_op = bbox(measured, baseline);

        auto full_bb = bb_op.value();

        point_list_t centers;
        for (size_t x_off = 0; x_off < 3; x_off++) {
            for (size_t y_off = 0; y_off < 3; y_off++) {
                auto bb = std::make_tuple(
                        std::get<0>(full_bb) + x_off * radius_size - 2 * radius_size,
                        std::get<1>(full_bb) + y_off * radius_size - 2 * radius_size,
                        std::get<2>(full_bb) + x_off * radius_size,
                        std::get<3>(full_bb) + y_off * radius_size);
                //Define a coarse grid.
                SparseGrid<wpt2_t> grid_red(bb, measured, radius_size * 3), grid_blue(bb, baseline, radius_size * 3);
                std::vector<uint64_t> keys;
                {
                    for (auto b = grid_red.begin(); b != grid_red.end(); b++) {
                        keys.emplace_back(b->first);
                    }
                    for (auto b = grid_blue.begin(); b != grid_blue.end(); b++) {
                        keys.emplace_back(b->first);
                    }
                    std::sort(keys.begin(), keys.end());
                    auto nend = std::unique(keys.begin(), keys.end());
                    keys.erase(nend, keys.end());
                }

                //Consider cells with points in them, so they must have $\eps n$ values
                for (auto k : keys) {

                    auto[lx, ly] = grid_red.get_lower_corner(k);
                    auto new_bb = std::make_tuple(
                            lx + radius_size,
                            ly + radius_size,
                            lx + 2 * radius_size,
                            ly + 2 * radius_size);
                    //std::cout << lx + radius_size << " " << ly + radius_size << " " << lx + 2* radius_size << " " << ly + 2 * radius_size << std::endl;
                    wpoint_list_t measured_local;
                    wpoint_list_t baseline_local;
                    auto[r_it, r_e_it] = grid_red(k);

                    double subgrid_red_total = 0.0;
                    for (auto b = r_it; b != r_e_it; b++) {
                        subgrid_red_total += b->second.get_weight();
                        measured_local.emplace_back(b->second);
                    }
                    auto[b_it, b_e_it] = grid_blue(k);
                    double subgrid_blue_total = 0.0;
                    for (auto b = b_it; b != b_e_it; b++) {
                        subgrid_blue_total += b->second.get_weight();
                        baseline_local.emplace_back(b->second);
                    }
                    double dense_grid_res = grid_res * (red_tot + blue_tot) / (subgrid_blue_total + subgrid_red_total);

                    max_kernel_slow_centers(
                            centers,
                            dense_grid_res,
                            new_bb);
                }
            }
        }
        return centers;
    }
}