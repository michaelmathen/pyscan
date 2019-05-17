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
        } while (status == GSL_CONTINUE && iter < 100);

        double p = gsl_vector_get(s->x, 0);
        double q = gsl_vector_get(s->x, 1);

        gsl_multimin_fdfminimizer_free (s);
        gsl_vector_free (x);
        double f_val = disc_f.lrt(p, q);

        return std::make_tuple(p, q, f_val);
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
        double max_v = -std::numeric_limits<double>::infinity();

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

                        disc_local->set_params(std::move(mr_j), std::move(br_j), radii);
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

    inline static bool colinear(const pt2_t &pt1, const pt2_t &pt2, const pt2_t &pt3) {
        double x1 = pt1(0), x2 = pt2(0), x3 = pt3(0);
        double y1 = pt1(1), y2 = pt2(1), y3 = pt3(1);

        return util::aeq(util::det2(x2 - x1, y2 - y1, x2 - x3, y2 - y3), 0.0);
    }

    inline static bool valid_pt(const pt2_t &p1, const pt2_t &p2, const pt2_t &p) {
        return !(p1.approx_eq(p) || p2.approx_eq(p) || colinear(p1, p2, p));
    }

    using annuli_t = std::vector<Disk>;

    inline static void get_annuli(
            const pt2_t& p1, const pt2_t& p2,
            const point_list_t& net,
            const std::vector<double>& res_ratio,
            double max_r,
            size_t curr_res,
            std::vector<annuli_t>& rings) {

        double orthoX = p2(1) - p1(1);
        double orthoY = p1(0) - p2(0);
        double cX = (p1(0) + p2(0)) / 2.0;
        double cY = (p1(1) + p2(1)) / 2.0;
        auto get_order = [orthoX, orthoY, cX, cY](const Disk &x) {
            auto origin = x.getOrigin();
            return orthoX * (origin(0) - cX) + orthoY * (origin(1) - cY) ;
        };

        rings.reserve(net.size());
        std::for_each(net.begin(), net.end(), [&](const pt2_t &p) {
            if (valid_pt(p1, p2, p)) {
                Disk tmp(p1, p2, p);
                double curr_res_ratio = tmp.getRadius() / res_ratio[curr_res];
                double outer_radius = curr_res_ratio * res_ratio.back();
                if (outer_radius <= max_r) {
                    annuli_t annulus;
                    for (auto r : res_ratio) {
                        auto center = tmp.getOrigin();
                        annulus.emplace_back(Disk(center(0), center(1), curr_res_ratio * r));
                    }
                    rings.emplace_back(annulus);
                }
            }
        });

        std::sort(rings.begin(), rings.end(), [&](const annuli_t &a, const annuli_t &b) {
            return get_order(a[0]) < get_order(b[0]);
        });
    }

    /*
     * We want to search a function that increases to a local max then decreases to a local min and then increases to a
     * new local max on the boundary of the function.
     *
     * So there are two places a point can exist inside of the disk. Either at that first value or at the new local min.
     *
     */
    template<typename It, typename Unary>
    It max_search(It begin, It end, Unary const& f)  {
//        auto exact_max = std::max_element(begin, end, [&f] (auto const& v1, auto const& v2){
//            return f(v1) < f(v2);
//        });
//
//        auto tmp_b = begin;
//        auto tmp_e = end;
        //If begin is monotonically increasing segment.
        while ((end - begin) >= 4) {
            auto mid = (end - begin) / 2 + begin;
            //In this case we are either the middle or to the left of the middle.
            if (f(*mid) < f(*(mid + 1))) {
                begin = mid + 1;
            } else if (f(*mid) > f(*(mid + 1))) {
                end = mid + 1;
            } else {
                break;
            }
        }
        auto fast_max = std::max_element(begin, end, [&f] (auto const& v1, auto const& v2){
            return f(v1) < f(v2);
        });
//        bool in_interval = false;
//        size_t interval_count = 0;
//        for (auto it = tmp_b; it != tmp_e; ++it) {
//            if (f(*it) <= 0) {
//                if (!in_interval) {
//                    interval_count += 1;
//                    in_interval = true;
//                }
//            } else {
//                in_interval = false;
//            }
//        }
//        if (interval_count == 2) {
//            std::cout << interval_count << " ";
//            for (auto it = tmp_b; it != tmp_e; ++it) {
//                std::cout << f(*it) << " ";
//            }
//            std::cout << std::endl;
//        }
//            for (auto it = tmp_b; it != tmp_e; ++it) {
//                std::cout << f(*it) << " ";
//            }
//            std::cout << std::endl;
//        }

        return fast_max;
    }

    /*
     * Sets begin to the first negative value and end to the first positive value or end of the range.
     * If it can't find a negative range then we set
     */
    template<typename It, typename Unary>
    std::tuple<It, It> find_range(It begin, It end, Unary const& f) {
        auto max_element = max_search(begin, end, f);
        if (max_element == end || f(*max_element) < 0) {
            return std::make_tuple(end, end);
        }

        auto interval_start = std::lower_bound(begin, max_element, 0.0, [&f] (auto const& v1, double v2) {
            return f(v1) < v2;
        });

        auto interval_end = std::upper_bound(max_element, end, 0.0, [&f] (double v1, auto const& v2) {
            return  v1 > f(v2);
        });
        return std::make_tuple(interval_start, interval_end);
    }

    //s points and roughly n disks with r annuli
    //idea 1. Find the mid disk first in log (n) time then scan to find left and right boundary for each disk.
    // takes O(n) time for each point.


    // Find a lower bound disk and an upper bound disk.
    std::vector<double> initialize_intervals(
            const pt2_t& p1,
            const pt2_t& p2,
            const wpoint_list_t& wpts,
            size_t curr_res,
            const std::vector<double>& annuli_scale,
            const std::vector<annuli_t>& annuli) {

        double orthoX = p2(1) - p1(1);
        double orthoY = p1(0) - p2(0);
        double scale = sqrt(orthoX * orthoX + orthoY * orthoY);

        double cX = (p1(0) + p2(0)) / 2.0;
        double cY = (p1(1) + p2(1)) / 2.0;

        pt2_t disk_seq_center(cX, cY, 1.0);
        double c_dist = p1.square_dist(p2) / 4.0;

        auto get_projected_pt = [orthoX, orthoY, scale, cX, cY](const pt2_t& x) {
            return pt2_t(orthoX * (x(0) - cX) / scale + cX, orthoY * (x(1) - cX) / scale + cY, 1.0);
        };

        std::vector<double> intervals(annuli.size(), 0.0);
        for (auto& wpt : wpts) {
            //Get the first disk that contains this point
            auto projected_pt = get_projected_pt(wpt);
            double thresh = 1 / annuli_scale[curr_res] * wpt.square_dist(projected_pt) - c_dist;
            // Two options.
            auto interval_dist = [&](const std::vector<Disk>& d_seq){
                auto& d1 = d_seq[curr_res];
                auto center = d1.getOrigin();
                return center.square_dist(disk_seq_center) - 1 / annuli_scale[curr_res] * projected_pt.square_dist(center)
                       - thresh;
            };
            auto [lb, ub] = find_range(annuli.begin(), annuli.end(), interval_dist);

            if (lb == ub) {
                //No disk in the sequence contains the pt.
                continue;
            }
            intervals[lb - annuli.begin()] += wpt.get_weight();
            if (ub != annuli.end()) {
                intervals[ub - annuli.begin()] -= wpt.get_weight();
            }
//            } else {
//                This is for if the scale is greater than 1, but I didn't want to handle this case.
//                auto interval_dist = [&](const std::vector<Disk>& d_seq){
//                    auto& d1 = d_seq[curr_res];
//                    auto center = d1.getOrigin();
//                    return thresh - center.square_dist(disk_seq_center) + 1 / annuli_scale[curr_res] * projected_pt.square_dist(center);
//                };
//                auto [lb, ub] = find_range(annuli.begin(), annuli.end(), interval_dist);
//                if (*lb == )
//                intervals[0] += wpt.get_weight();
//                intervals[lb - annuli.begin()] -= wpt.get_weight();
//                intervals[ub - annuli.begin()] += wpt.get_weight();
//            }
        }
        return intervals;
    }

    inline static std::tuple<Disk, double> max_annuli_multi(
            const pt2_t &p1, const pt2_t &p2,
            const point_list_t &net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            const std::vector<double>& annuli_res,
            double max_r,
            KDisc& disc) {

        Disk max_disk;
        double max_v = 0;
        if (p1.approx_eq(p2)) {
            return std::make_tuple(max_disk, 0.0);
        }

        std::vector<annuli_t> rings;
        // Check to see how slow this is. We might be losing a lot of time, by sorting these big objects.
        // Sort a set of ordered objects first and then reorder afterwards.
        get_annuli(p1, p2, net, annuli_res, max_r, annuli_res.size() - 1, rings);

        if (rings.empty()) {
            return std::make_tuple(max_disk, 0.0);
        }

        std::vector<std::vector<double>> red_resolution_deltas(annuli_res.size());
        std::vector<std::vector<double>> blue_resolution_deltas(annuli_res.size());
        for (size_t i = 0; i < annuli_res.size(); i++) {
            red_resolution_deltas[i] = initialize_intervals(p1, p2, red, i, annuli_res, rings);
            blue_resolution_deltas[i] = initialize_intervals(p1, p2, blue, i, annuli_res, rings);
        }

        std::vector<double> curr_radii(annuli_res.size());
        std::vector<double> curr_mr(annuli_res.size(), 0.0);
        std::vector<double> curr_br(annuli_res.size(), 0.0);

        double p = .6, q = .5;
        for (size_t i = 0; i < rings.size(); i++) {
            //Initialize the current set of radii.
            for (size_t j = 0; j < annuli_res.size(); j++){
                curr_radii[j] = rings[i][j].getRadius();
                curr_mr[j] += red_resolution_deltas[j][i];
                curr_br[j] += blue_resolution_deltas[j][i];
            }
            disc.set_params(curr_mr, curr_br, curr_radii);
            double fval;
            std::tie(p, q, fval) = find_pq_poi(p, q, disc);
            if (std::isnan(p) || std::isnan(q) || 0 >= p || p >= 1 || 0 >= q || q >= 1) {
                p = .6;
                q = .5;
            }
            if (max_v < fval) {
                max_v = fval;
                //Set the biggest disk
                max_disk = rings[i].back();
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
        double max_v = -std::numeric_limits<double>::infinity();
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
        double max_stat = -std::numeric_limits<double>::infinity();
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

            if (net_chunk.size() >= 2) {
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

    std::tuple<Disk, double> max_annuli_scale_multi(
            const point_list_t &point_net,
            const std::vector<wpt2_t> &red,
            const std::vector<wpt2_t> &blue,
            std::vector<double> res_scales,
            double max_radii,
            const KDisc& disc) {

        Disk cur_max;
        double max_stat = 0.0;
        if (res_scales.empty()) {
            return std::make_tuple(cur_max, max_stat);
        }
        if (point_net.empty()) {
            return std::make_tuple(Disk(), 0.0);
        }
        auto bb_op = bbox(point_net, red, blue);
        if (!bb_op.has_value()) {
            return std::make_tuple(cur_max, max_stat);
        }
        std::sort(res_scales.begin(), res_scales.end());
        for (size_t i = 0; i < res_scales.size(); i++) {
            res_scales[i] = res_scales[i] / res_scales.back();
        }

        auto disc_local = disc.get_copy();
        auto bb = bb_op.value();
        SparseGrid<pt2_t> grid_net(bb, point_net, max_radii);
        auto grid_r = grid_net.get_grid_size();
        SparseGrid<wpt2_t> grid_red(bb, red, max_radii), grid_blue(bb, blue, max_radii);


        for (auto center_cell = grid_net.begin(); center_cell != grid_net.end();) {
            std::vector<pt2_t> net_chunk;
            std::vector<wpt2_t> red_chunk;
            std::vector<wpt2_t> blue_chunk;
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
                    for (auto &pt2: net_chunk) {
                        if (pt1->second.approx_eq(pt2)) continue;

                        auto[local_max_disk, local_max_stat] =
                        max_annuli_multi(pt1->second, pt2, net_chunk, red_chunk, blue_chunk,
                                res_scales, max_radii, *(disc_local.get()));
                        if (local_max_stat > max_stat) {
                            cur_max = local_max_disk;
                            max_stat = local_max_stat;
                        }
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