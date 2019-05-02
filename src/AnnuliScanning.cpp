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

#include "AnnuliScanning.hpp"
#include "SparseGrid.hpp"
#include "Utilities.hpp"


namespace pyscan {


    kernel_func_t gauss_kernel(double deviation){
        return [deviation] (double dist) {
            return 1 / sqrt(2 * M_PI) * exp(- dist * dist / (2 * deviation * deviation));
        };
    }


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
        auto fp = static_cast<KDisc*>(params);

        double p = gsl_vector_get(v, 0);
        double q = gsl_vector_get(v, 1);

        return fp->get_function()(p, q);
    }

    /* The gradient of f, df = (df/dx, df/dy). */
    void
    unwrap_df (const gsl_vector *v, void *params,
           gsl_vector *df)
    {
        auto fp = static_cast<KDisc*>(params);

        double p = gsl_vector_get(v, 0);
        double q = gsl_vector_get(v, 1);

        auto [dfp, dfq] = fp->get_differential()(p, q);
        gsl_vector_set(df, 0, dfp);
        gsl_vector_set(df, 1, dfq);
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
            const discrepancy_kfunc_t& disc_f) {
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
                status = gsl_multimin_test_gradient(s->gradient, 1e-5);
            }
            p_o = p;
            q_o = q;
            //std::cout << "p = " << gsl_vector_get(s->x, 0) << " q = " << gsl_vector_get(s->x, 1) << "f = " << s->f <<std::endl;
        } while (status == GSL_CONTINUE && iter < 100);

        std::cout << std::endl;
        double p = gsl_vector_get(s->x, 0);
        double q = gsl_vector_get(s->x, 1);

        gsl_multimin_fdfminimizer_free (s);
        gsl_vector_free (x);
        return std::make_tuple(p, q, s->f);
    }

    std::vector<double> propagate_annuli(const wpoint_list_t& pts, const pt2_t & center, const std::vector<double>& radii) {
        size_t i = 0;
        std::vector<double> valr_j(radii.size(), 0.0);
        for (auto& p : pts) {
            auto curr_dist = center.dist(p);
            while (i < radii.size() && curr_dist > radii[i]) {
                i++;
            }
            if (i >= radii.size()) {
                break;
            }
            if (center.dist(p) < radii[i]) {
                valr_j[i] += p.get_weight();
            }
        }
        return valr_j;
    }

    std::tuple<Disk, double> max_annuli(const point_list_t &pts,
                                        wpoint_list_t mpts,
                                        wpoint_list_t bpts,
                                        const std::vector<double> &radii,
                                        const KDisc& disc) {


        double p_init = .6, q_init = .5;
        Disk max_disk;
        double max_v = 0;

        if (pts.size() < 2 || radii.empty()) {
            return std::make_tuple(max_disk, max_v);
        }

        auto disc_local = disc.get_copy();
        for (auto r_it = radii.begin(); r_it != radii.end(); ++r_it) {

            for (auto nit1 = pts.begin(); nit1 != pts.end() - 1; ++nit1) {
                for (auto nit2 = nit1 + 1; nit2 != pts.end(); ++nit2) {
                    // Check to make sure the points are close enough to be on the boundary of some disk of radii
                    // r
                    if (nit1->square_dist(*nit2) < 4 * (*r_it) * (*r_it)) {
                        //Now we can get the origin point and compute all the annuli
                        Disk test_disk(*nit1, *nit2, *r_it);
                        auto center = test_disk.getOrigin();

                        //Now we have the annuli
                        std::sort(mpts.begin(), mpts.end(), [&](const pt2_t& p1, const pt2_t& p2) {
                            return center.square_dist(p1) < center.square_dist(p2);
                        });

                        std::sort(bpts.begin(), bpts.end(), [&](const pt2_t& p1, const pt2_t& p2) {
                            return center.square_dist(p1) < center.square_dist(p2);
                        });
                        // Accumulate to an annulus
                        auto mr_j = propagate_annuli(mpts, center, radii);
                        auto br_j = propagate_annuli(bpts, center, radii);

                        disc_local->set_params(mr_j, br_j, radii);
//                        auto fval = disc_local->get_function()(p_init, q_init);
//                        std::cout << fval << std::endl;
//                        std::cout << mr_j << " " << br_j << std::endl;

                        auto [p, q, fval] = find_pq_poi(p_init, q_init, *(disc_local.get()));
                        if (max_v < fval) {
                            max_v = fval;
                            max_disk = Disk(center(0), center(1), *r_it);
                        }
                    }
                }
            }
        }
        return std::make_tuple(max_disk, max_v);
    }

    std::tuple<Disk, double> max_annuli_restricted(
            pt2_t const& pt,
            const point_list_t &pts,
            wpoint_list_t mpts,
            wpoint_list_t bpts,
            const std::vector<double> &radii,
            const KDisc &disc) {

        double p_init = .6, q_init = .5;
        Disk max_disk;
        double max_v = 0;
        auto disc_local = disc.get_copy();
        for (auto r_it = radii.begin(); r_it != radii.end(); ++r_it) {

            for (auto nit1 = pts.begin(); nit1 != pts.end() - 1; ++nit1) {
                // Check to make sure the points are close enough to be on the boundary of some disk of radii
                // r
                if (pt.square_dist(*nit1) < 4 * (*r_it) * (*r_it)) {
                    //Now we can get the origin point and compute all the annuli
                    Disk test_disk(*nit1, pt, *r_it);
                    auto center = test_disk.getOrigin();

                    //Now we have the annuli
                    std::sort(mpts.begin(), mpts.end(), [&](const pt2_t &p1, const pt2_t &p2) {
                        return center.square_dist(p1) < center.square_dist(p2);
                    });

                    std::sort(bpts.begin(), bpts.end(), [&](const pt2_t &p1, const pt2_t &p2) {
                        return center.square_dist(p1) < center.square_dist(p2);
                    });
                    auto mr_j = propagate_annuli(mpts, center, radii);
                    auto br_j = propagate_annuli(bpts, center, radii);
                    disc_local->set_params(mr_j, br_j, radii);
                    auto [p, q, fval] = find_pq_poi(p_init, q_init, *(disc_local.get()));
                    if (max_v < fval) {
                        max_v = fval;
                        max_disk = Disk(center(0), center(1), *r_it);
                    }
                }
            }
        }

        return std::make_tuple(max_disk, max_v);
    }


    std::tuple<Disk, double> max_annuli_scale(
            const point_list_t &point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            const std::vector<double>& annuli_res,
            const KDisc& f) {

        Disk cur_max;
        double max_stat = 0.0;
        if (point_net.empty() || annuli_res.empty()) {
            return std::make_tuple(Disk(), 0.0);
        }
        auto bb_op = bbox(point_net, red, blue);
        if (!bb_op.has_value()) {
            return std::make_tuple(cur_max, max_stat);
        }
        auto bb = bb_op.value();
        SparseGrid<pt2_t> grid_net(bb, point_net, annuli_res.back());
        auto grid_r = grid_net.get_grid_size();
        SparseGrid<wpt2_t> grid_red(bb, red, annuli_res.back()), grid_blue(bb, blue, annuli_res.back());


        for (auto center_cell = grid_net.begin(); center_cell != grid_net.end();) {
            std::vector<pt2_t> net_chunk;
            wpoint_list_t red_chunk;
            wpoint_list_t blue_chunk;
            net_chunk.clear();
            red_chunk.clear();
            blue_chunk.clear();
            size_t i, j;
            std::tie(i, j) = grid_net.get_cell(center_cell->second);
            size_t start_k = i < 2 ? 0 : i - 2;
            size_t start_l = j < 2 ? 0 : j - 2;
            size_t end_k = i + 2 < grid_r ? i + 2 : grid_r;
            size_t end_l = j + 2 < grid_r ? j + 2 : grid_r;
            auto range = grid_net(i, j);

            for (size_t k = start_k; k <= end_k; ++k) {
                for (size_t l = start_l; l <= end_l; ++l) {
                    auto net_range = grid_net(k, l);
                    for (auto it = net_range.first; it != net_range.second; ++it) {
                        net_chunk.emplace_back(it->second);
                    }

                    auto red_range = grid_red(k, l);
                    for (auto it = red_range.first; it != red_range.second; ++it)
                        red_chunk.emplace_back(it->second);


                    auto blue_range = grid_blue(k, l);
                    for (auto it = blue_range.first; it != blue_range.second; ++it)
                        blue_chunk.emplace_back(it->second);
                }
            }

            if (net_chunk.size() >= 3) {
                for (auto pt1 = range.first; pt1 != range.second; ++pt1) {
                    auto [local_max_disk, local_max_stat] =
                    max_annuli_restricted(pt1->second, net_chunk, red_chunk, blue_chunk,
                                        annuli_res,f);
                    if (local_max_stat > max_stat) {
                        cur_max = local_max_disk;
                        max_stat = local_max_stat;
                    }
                }
            }


            auto last = center_cell->first;
            do {
                ++center_cell;
            } while (center_cell != grid_net.end() && center_cell->first == last);
        }

        return std::make_tuple(cur_max, max_stat);
    }

}