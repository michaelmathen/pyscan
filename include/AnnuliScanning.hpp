/*
 * Created by Michael Matheny on 4/25/19.
 * at the University of Utah
 * email: michaelmathen@gmail.com
 * website: https://mmath.dev/
 */

#ifndef PYSCAN_ANNULISCANNING_HPP
#define PYSCAN_ANNULISCANNING_HPP

#include "Disk.hpp"
#include "Point.hpp"

namespace pyscan {

    using disk_list_t = std::vector<Disk>;
    using kernel_func_t = std::function<double(double, double)>;
    using pq_disc_t = std::function<double(double, double)>;
    using pqdf_disc_t = std::function<std::tuple<double, double>(double, double)>;


    //generate a guassian kernel
//    kernel_func_t rectangle_kernel(double deviation);
//
//    kernel_func_t triangle_kernel(double deviation);
//
//    kernel_func_t epanechnikov_kernel(double deviation);

    struct KDisc {

        pq_disc_t f;
        pqdf_disc_t df;
        kernel_func_t kern;
        double bandwidth_ratio;
        std::vector<double> m_annuli;
        std::vector<double> b_annuli;
        std::vector<double> radii;

        explicit KDisc(kernel_func_t kernelFunc, double br) : kern(std::move(kernelFunc)), bandwidth_ratio(br) {}

        virtual void set_params(
                std::vector<double> m_annuli,
                std::vector<double> b_annuli,
                std::vector<double> radii) = 0;

        virtual pq_disc_t& get_function() {
            return f;
        }

        virtual pqdf_disc_t& get_differential() {
            return df;
        }
        virtual std::shared_ptr<KDisc> get_copy() const = 0;

        virtual double lrt(double p, double q) const = 0;
    };

    inline double bernoulli_f(
            double p, double q,
            const std::vector<double>& m_annuli,
            const std::vector<double>& b_annuli,
            const std::vector<double>& radii,
            double bandwidth,
            const kernel_func_t& kernel) {
        double val = 0;
        double prev_dist = 0;
        for (size_t i = 0; i < radii.size(); i++) {
            double fr = kernel((prev_dist + radii[i]) / 2.0, bandwidth);
            double gr = fr * p + (1 - fr) * q;
            val += m_annuli[i] * log(gr) + b_annuli[i] * log(1 - gr);
            prev_dist = radii[i];
        }
        return -val;
    }

    inline std::tuple<double, double> bernoulli_df(
            double p, double q,
            const std::vector<double>& m_annuli,
            const std::vector<double>& b_annuli,
            const std::vector<double>& radii,
            double bandwidth,
            const kernel_func_t& kernel) {

        double dp = 0.0, dq = 0.0;
        double prev_dist = 0;

        for (size_t i = 0; i < radii.size(); i++) {
            double fr = kernel((prev_dist + radii[i]) / 2.0, bandwidth);

            dp += -m_annuli[i] * fr / ((p - q) * fr + q);
            dp += b_annuli[i] * fr / ((1 - q) + (q - p) * fr);

            dq += -m_annuli[i] * (1 - fr) / ((p - q) * fr + q);
            dq += b_annuli[i] * (1 - fr) / ((1 - q) + (q - p) * fr);
            prev_dist = radii[i];
        }
        return std::make_tuple(dp, dq);
    }

    struct Bernouli_kf : public KDisc {

        Bernouli_kf(kernel_func_t kernel, double bandwidth_ratio) : KDisc(std::move(kernel), bandwidth_ratio) {}
        void set_params(
                std::vector<double> m_annuli,
                std::vector<double> b_annuli,
                std::vector<double> radii) final {
            double bandwidth = radii.back() * bandwidth_ratio;
            f = [&, bandwidth](double p, double q) {
                return bernoulli_f(p, q, m_annuli, b_annuli, radii, bandwidth, kern);
            };
            df = [&, bandwidth](double p, double q) {
                return bernoulli_df(p, q, m_annuli, b_annuli, radii, bandwidth, kern);
            };
            this->m_annuli = m_annuli;
            this->b_annuli = b_annuli;
            this->radii = radii;

        }

        std::shared_ptr<KDisc> get_copy() const final {
            auto ptr = std::shared_ptr<KDisc>(new Bernouli_kf(this->kern, bandwidth_ratio));
            ptr->kern = this->kern;
            return ptr;
        }

        double lrt(double p, double q) const final {
            double prev_dist = 0;
            double lrt_val = 0;
            double mr = 0;
            double br = 0;
            for (size_t i = 0; i < radii.size(); i++) {
                mr += m_annuli[i];
                br += b_annuli[i];
            }
            double scale = mr / (mr + br);
            double bandwidth = bandwidth_ratio * radii.back();
            for (size_t i = 0; i < radii.size(); i++) {
                double fr = kern((prev_dist + radii[i]) / 2.0, bandwidth);
                double gr = fr * p + (1 - fr) * q;
                lrt_val += m_annuli[i] * log(gr / scale) + b_annuli[i] * log((1 - gr) / (1 - scale));
                prev_dist = radii[i];
            }
            return lrt_val;
        }
    };


    using discrepancy_kfunc_t = KDisc;


    std::tuple<Disk, double> max_annuli(const point_list_t &pts,
                                        wpoint_list_t mpts,
                                        wpoint_list_t bpts,
                                        const std::vector<double> &radii,
                                        const KDisc& disc);


    std::tuple<Disk, double> max_annuli_restricted(
            pt2_t const& pt,
            const point_list_t &pts,
            wpoint_list_t mpts,
            wpoint_list_t bpts,
            const std::vector<double> &radii,
            const KDisc &disc);

    std::tuple<Disk, double> max_annuli_scale(
            const point_list_t &point_net,
            const wpoint_list_t &red,
            const wpoint_list_t &blue,
            const std::vector<double>& annuli_res,
            const KDisc& f);

    std::tuple<Disk, double> max_annuli_scale_multi(
            const point_list_t &point_net,
            const std::vector<wpt2_t> &red,
            const std::vector<wpt2_t> &blue,
            std::vector<double> res_scales,
            double max_radii,
            const KDisc& disc);
}
#endif //PYSCAN_ANNULISCANNING_HPP
