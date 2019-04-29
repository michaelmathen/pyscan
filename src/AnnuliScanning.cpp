/*
 * Created by Michael Matheny on 4/25/19.
 * at the University of Utah
 * email: michaelmathen@gmail.com
 * website: https://mmath.dev/
 */

#include "AnnuliScanning.hpp"
#include "Utilities.hpp"

namespace pyscan {
    using disk_list_t = std::vector<Disk>;

    class WDisk : public Disk {
        double weight;
    public:
        WDisk(double w_, double x_, double y_, double r_)
                : Disk(x, y, r_), weight(w_) {}
        WDisk()
                : Disk(), weight(0.0) {}
        WDisk(double w_, const pt2_t &pt1, const pt2_t &pt2, const pt2_t &pt3) : Disk(pt1, pt2, pt3), weight(w_) {}

        WDisk(double w_, const pt2_t &p1, const pt2_t &p2, double r_) : Disk(p1, p2, r_), weight(w_) {}

        double get_weight() const {
            return weight;
        }
    };
    using wdisk_list_t = std::vector<WDisk>;


    void incr_it(wpoint_it_t b, wpoint_it_t e, wpoint_it_t& w_it){
        size_t i = ((w_it - b) + 1) % (e - b);
        w_it = b + i;
    }

    static inline double order_f(const Point<>& p1, const Disk& d1) {
        return p1.direction(d1.getOrigin())[0];
    }


    std::tuple<wpoint_it_t, wpoint_it_t> increment_window(
            wpoint_it_t begin,
            wpoint_it_t end,
            wpoint_it_t w_begin,
            wpoint_it_t w_end,
            const pt2_t& center,
            const Disk& d_new) {
        auto n_end = w_end;
        //The new end is either inside the window or before the window.
        auto ahead = n_end;
        incr_it(begin, end, ahead);
        for (; ahead != w_end; incr_it(begin, end, ahead)) {
            Disk disk(*ahead, center, d_new.getRadius());
            if (order_f(center, d_new) < order_f(center, disk)) {
                break;
            }
            n_end = ahead;
        }

        //The new begining is either inside the window or before the window.
        //If it is inside then we are done.
        auto n_begin = w_begin;
        if (!d_new.contains(*n_begin)) {
            incr_it(begin, end, n_begin);
            for (; n_begin != n_end; incr_it(begin, end, w_begin)) {
                if (d_new.contains(*n_begin)) {
                    break;
                }
            }
        }
        return std::make_tuple(n_begin, n_end);
    }

    std::tuple<wpoint_it_t, wpoint_it_t> initial_window(
            wpoint_it_t begin,
            wpoint_it_t end,
            const pt2_t& center,
            const Disk& d_new){
        wpoint_it_t n_end;
        auto ahead = begin;
        incr_it(begin, end, ahead);
        for (; ahead != begin; incr_it(begin, end, ahead)) {
            Disk disk(*ahead, center, d_new.getRadius());
            if (order_f(center, d_new) < order_f(center, disk)) {
                break;
            }
            n_end = ahead;
        }

        for (; ahead != n_end; incr_it(begin, end, ahead)) {
            if (d_new.contains(*ahead)) {
                break;
            }
        }
        return std::make_tuple(ahead, n_end);
    }

    double sum_interval(wpoint_it_t b, wpoint_it_t e) {
        double sum = 0;
        for (; b != e; b++) {
            sum += b->get_weight();
        }
        return sum;
    }

    std::tuple<Disk, double> max_annuli(const point_list_t &pts,
                                        wpoint_list_t mpts,
                                        wpoint_list_t bpts,
                                        const std::vector<double> &radii,
                                        const discrepancy_func_t &func) {
        for (auto r_it = radii.begin(); r_it != radii.end(); ++r_it) {

            for (auto nit1 = pts.begin(); nit1 != pts.end() - 1; ++nit1) {
                for (auto nit2 = nit1 + 1; nit2 != pts.end(); ++nit2) {
                    // Check to make sure the points are close enough to be on the boundary of some disk of radii
                    // r
                    if (nit1->square_dist(*nit2) < 4 * (*r_it) * (*r_it)) {
                        net_disks.emplace_back(*nit1, *nit2, *r_it);
                    }
                }
            }
        }
    }

//    std::tuple<Disk, double> max_annuli(const point_list_t &pts,
//                                        wpoint_list_t mpts,
//                                        wpoint_list_t bpts,
//                                        const std::vector<double> &radii,
//                                        const discrepancy_func_t &func) {
//
//        for (auto r_it = radii.begin(); r_it != radii.end(); ++r_it) {
//            for (auto nit1 = pts.begin(); nit1 != pts.end() - 1; ++nit1) {
//
//                disk_list_t net_disks;
//                for (auto nit2 = nit1 + 1; nit2 != pts.end(); ++nit2) {
//                    // Check to make sure the points are close enough to be on the boundary of some disk of radii
//                    // r
//                    if (nit1->square_dist(*nit2) < 4 * (*r_it) * (*r_it)) {
//                        net_disks.emplace_back(*nit1, *nit2, *r_it);
//                    }
//                }
//                if (net_disks.size() < 1) {
//                    break;
//                }
//
//                //Sort by the orientation with the initial point.
//                std::sort(net_disks.begin(), net_disks.end(), [&](const Disk& d1, const Disk& d2) {
//                    return order_f(*nit1, d1) < order_f(*nit1, d2);
//                });
//
//                auto d0 = net_disks[0];
//                auto order_and_partition = [&](wpoint_list_t & wpts) {
//                    std::vector<double> ordering;
//                    auto nend = std::partition(wpts.begin(), wpts.end(), [&](const pt2_t& p1) {
//                        return nit1->square_dist(p1) < 4 * (*r_it) * (*r_it);
//                    });
//                    std::sort(wpts.begin(), nend, [&](const pt2_t& p1, const pt2_t& p2) {
//                        Disk d1(*nit1, p1, *r_it);
//                        Disk d2(*nit1, p2, *r_it);
//                        return order_f(*nit1, d1) < order_f(*nit1, d2);
//                    });
//                    return nend;
//                };
//
//                auto m_end = order_and_partition(mpts);
//                auto b_end = order_and_partition(bpts);
//                auto [mb, me] = initial_window(mpts.begin(), m_end, *nit1, d0);
//                auto [bb, be] = initial_window(bpts.begin(), b_end, *nit1, d0);
//                double m_count = sum_interval(mb, me);
//                double b_count = sum_interval(bb, be);
//                for (auto& disk : net_disks) {
//                    auto [nmb, nme] = increment_window(mpts.begin(), m_end, mb, me, *nit1, disk);
//                    auto [nbb, nbe] = increment_window(bpts.begin(), b_end, bb, be, *nit1, disk);
//
//                    m_count += sum_interval(me, nme) - sum_interval(mb, nmb);
//                    b_count += sum_interval(be, nbe) - sum_interval(bb, nbb);
//                }
//
//
//            }
//        }
//    }
}