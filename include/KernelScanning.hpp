/*
 * Created by Michael Matheny on 4/25/19.
 * at the University of Utah
 * email: michaelmathen@gmail.com
 * website: https://mmath.dev/
 */

#ifndef PYSCAN_KERNELSCANNING_HPP
#define PYSCAN_KERNELSCANNING_HPP

#include "Disk.hpp"
#include "Point.hpp"

namespace pyscan {

//    using disk_list_t = std::vector<Disk>;
    using kernel_func_t = std::function<double(double, double)>;
//    using pq_disc_t = std::function<double(double, double)>;
//    using pqdf_disc_t = std::function<std::tuple<double, double>(double, double)>;


    //generate a guassian kernel
//    kernel_func_t rectangle_kernel(double deviation);
//
//    kernel_func_t triangle_kernel(double deviation);
//
//    kernel_func_t epanechnikov_kernel(double deviation);

//    struct KDisc {
//
//        pq_disc_t f;
//        pqdf_disc_t df;
//        kernel_func_t kern;
//        double bandwidth_ratio;
//        std::vector<double> m_annuli;
//        std::vector<double> b_annuli;
//        std::vector<double> radii;
//
//        explicit KDisc(kernel_func_t kernelFunc, double br) : kern(std::move(kernelFunc)), bandwidth_ratio(br) {}
//
//        virtual void set_params(
//                std::vector<double> m_annuli,
//                std::vector<double> b_annuli,
//                std::vector<double> radii) = 0;
//
//        virtual pq_disc_t& get_function() {
//            return f;
//        }
//
//        virtual pqdf_disc_t& get_differential() {
//            return df;
//        }
//        virtual std::shared_ptr<KDisc> get_copy() const = 0;
//
//        virtual double lrt(double p, double q) const = 0;
//    };
//
//    inline double bernoulli_f(
//            double p, double q,
//            const std::vector<double>& m_annuli,
//            const std::vector<double>& b_annuli,
//            const std::vector<double>& radii,
//            double bandwidth,
//            const kernel_func_t& kernel) {
//        double val = 0;
//        double prev_dist = 0;
//        for (size_t i = 0; i < radii.size(); i++) {
//            double fr = kernel((prev_dist + radii[i]) / 2.0, bandwidth);
//            double gr = fr * p + (1 - fr) * q;
//            val += m_annuli[i] * log(gr) + b_annuli[i] * log(1 - gr);
//            prev_dist = radii[i];
//        }
//        return -val;
//    }
//
//    inline std::tuple<double, double> bernoulli_df(
//            double p, double q,
//            const std::vector<double>& m_annuli,
//            const std::vector<double>& b_annuli,
//            const std::vector<double>& radii,
//            double bandwidth,
//            const kernel_func_t& kernel) {
//
//        double dp = 0.0, dq = 0.0;
//        double prev_dist = 0;
//
//        for (size_t i = 0; i < radii.size(); i++) {
//            double fr = kernel((prev_dist + radii[i]) / 2.0, bandwidth);
//
//            dp += -m_annuli[i] * fr / ((p - q) * fr + q);
//            dp += b_annuli[i] * fr / ((1 - q) + (q - p) * fr);
//
//            dq += -m_annuli[i] * (1 - fr) / ((p - q) * fr + q);
//            dq += b_annuli[i] * (1 - fr) / ((1 - q) + (q - p) * fr);
//            prev_dist = radii[i];
//        }
//        return std::make_tuple(dp, dq);
//    }
//
//    struct Bernouli_kf : public KDisc {
//
//        Bernouli_kf(kernel_func_t kernel, double bandwidth_ratio) : KDisc(std::move(kernel), bandwidth_ratio) {}
//        void set_params(
//                std::vector<double> m_annuli,
//                std::vector<double> b_annuli,
//                std::vector<double> radii) final {
//            double bandwidth = radii.back() * bandwidth_ratio;
//            f = [&, bandwidth](double p, double q) {
//                return bernoulli_f(p, q, m_annuli, b_annuli, radii, bandwidth, kern);
//            };
//            df = [&, bandwidth](double p, double q) {
//                return bernoulli_df(p, q, m_annuli, b_annuli, radii, bandwidth, kern);
//            };
//            this->m_annuli = m_annuli;
//            this->b_annuli = b_annuli;
//            this->radii = radii;
//
//        }
//
//        std::shared_ptr<KDisc> get_copy() const final {
//            auto ptr = std::shared_ptr<KDisc>(new Bernouli_kf(this->kern, bandwidth_ratio));
//            ptr->kern = this->kern;
//            return ptr;
//        }
//
//        double lrt(double p, double q) const final {
//            double prev_dist = 0;
//            double lrt_val = 0;
//            double mr = 0;
//            double br = 0;
//            for (size_t i = 0; i < radii.size(); i++) {
//                mr += m_annuli[i];
//                br += b_annuli[i];
//            }
//            double scale = mr / (mr + br);
//            double bandwidth = bandwidth_ratio * radii.back();
//            for (size_t i = 0; i < radii.size(); i++) {
//                double fr = kern((prev_dist + radii[i]) / 2.0, bandwidth);
//                double gr = fr * p + (1 - fr) * q;
//                lrt_val += m_annuli[i] * log(gr / scale) + b_annuli[i] * log((1 - gr) / (1 - scale));
//                prev_dist = radii[i];
//            }
//            return lrt_val;
//        }
//    };
//
//
//    using discrepancy_kfunc_t = KDisc;


//    std::tuple<Disk, double> max_annuli(const point_list_t &pts,
//                                        wpoint_list_t mpts,
//                                        wpoint_list_t bpts,
//                                        const std::vector<double> &radii,
//                                        const KDisc& disc);
//
//
//    std::tuple<Disk, double> max_annuli_restricted(
//            pt2_t const& pt,
//            const point_list_t &pts,
//            wpoint_list_t mpts,
//            wpoint_list_t bpts,
//            const std::vector<double> &radii,
//            const KDisc &disc);
//
//    std::tuple<Disk, double> max_annuli_scale(
//            const point_list_t &point_net,
//            const wpoint_list_t &red,
//            const wpoint_list_t &blue,
//            const std::vector<double>& annuli_res,
//            const KDisc& f);
//
//    std::tuple<Disk, double> max_annuli_scale_multi(
//            const point_list_t &point_net,
//            const std::vector<wpt2_t> &red,
//            const std::vector<wpt2_t> &blue,
//            std::vector<double> res_scales,
//            double max_radii,
//            const KDisc& disc);



    class Bernoulli_Disk {
    public:
        Bernoulli_Disk(double m_tot,
                        double b_tot,
                        double b,
                        const kernel_func_t& k) : m_total(m_tot), b_total(b_tot),
                bandwidth(b), kernel(k) {}


        double alternative_hyp(double p, double q) const {
            double m_total = this->m_total;
            double b_total = this->b_total;
            double val = 0;
            for (size_t i = 0; i < mr.size(); i++) {
                double fr = kernel(m_radii[i], bandwidth);
                //std::cout << fr << " " << m_radii[i] << std::endl;
                double gr = fr * p + (1 - fr) * q;
                val += mr[i] * log(gr);
                m_total -= mr[i];
            }
            for (size_t i = 0; i < br.size(); i++) {
                double fr = kernel(b_radii[i], bandwidth);
                double gr = fr * p + (1 - fr) * q;
                val += br[i] * log(1 - gr);
                b_total -= br[i];
            }
            val += m_total * log(q) + b_total * log(1 - q);

            return val;
        }

        double null_hyp() const {
            double scale = m_total / (m_total + b_total);
            return m_total * log(scale) + b_total * log(1 - scale);
        }


        inline std::tuple<double, double> mass_conserve(double p, double q) const {
            double m_total_tmp = this->m_total;
            double b_total_tmp = this->b_total;
            double m_total = this->m_total;
            double b_total = this->b_total;


            double mass_sum1 = 0.0;
            double mass_sum2 = 0.0;
            for (size_t i = 0; i < mr.size(); i++) {
                double fr = kernel(m_radii[i], bandwidth);
                double gr = fr * p + (1 - fr) * q;
                mass_sum1 += mr[i] * 1 / gr;
                m_total -= mr[i];
            }

            for (size_t i = 0; i < br.size(); i++) {
                double fr = kernel(b_radii[i], bandwidth);
                double gr = fr * p + (1 - fr) * q;
                mass_sum2 += br[i] * 1/(1 - gr);
                b_total -= br[i];
            }
            mass_sum1 += m_total / q;
            mass_sum2 += b_total / (1 - q);
            return std::make_tuple(mass_sum1 - m_total_tmp - b_total_tmp, mass_sum2 - b_total_tmp - m_total_tmp);
        }
        inline std::tuple<double, double> diff(double p, double q) const {
            double m_total = this->m_total;
            double b_total = this->b_total;
            double dp = 0.0;
            double dq = 0.0;
            for (size_t i = 0; i < mr.size(); i++) {
                double fr = kernel(m_radii[i], bandwidth);
                double gr = fr * p + (1 - fr) * q;
                dp += mr[i] * fr / gr;
                dq += mr[i] * (1 - fr) / gr;

                m_total -= mr[i];
            }

            for (size_t i = 0; i < br.size(); i++) {
                double fr = kernel(b_radii[i], bandwidth);
                double gr = fr * p + (1 - fr) * q;
                dp += -br[i] * fr / (1 - gr);
                dq += -br[i] * (1 - fr) / (1 - gr);
                b_total -= br[i];
            }
            dq += m_total / q - b_total / (1 - q);
            return std::make_tuple(dp, dq);
        }

        inline double lrt(double p, double q) const {
            return alternative_hyp(p, q)- null_hyp() ;
        }

        void set_weights(std::vector<double> m,
                       std::vector<double> b) {
            mr = std::move(m);
            br = std::move(b);
        }

        void set_radii(std::vector<double> mr_temp,
                         std::vector<double> br_temp) {
            m_radii = std::move(mr_temp);
            b_radii = std::move(br_temp);
        }


        std::vector<double> mr;
        std::vector<double> br;
        std::vector<double> m_radii;
        std::vector<double> b_radii;
        double m_total;
        double b_total;
        double bandwidth;
        kernel_func_t kernel;
    };

    std::tuple<double, double, double> measure_kernel(
            const pt2_t& center,
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double bandwidth);

    std::tuple<Disk, double> max_kernel_slow(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double disk_r,
            double bandwidth);

    point_list_t kernel_centers_approximate(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth);

    std::tuple<Disk, double> max_kernel_slow2(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth);

    std::tuple<Disk, double> max_kernel_adaptive(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth);

    std::tuple<Disk, double> max_kernel_prune_far(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth);

    std::tuple<Disk, double> max_kernel(
            const wpoint_list_t &measured,
            const wpoint_list_t &baseline,
            double grid_res,
            double radius_size,
            double bandwidth);

}
#endif //PYSCAN_KERNELSCANNING_HPP
