//
// Created by mmath on 9/25/17.
//
#include <algorithm>
#include <unordered_map>
#include <ostream>
#include <iterator>
#include <set>
#include <deque>
#include <unordered_set>
#include <cassert>

#include "Range.hpp"
#include "SparseGrid.hpp"
#include "Statistics.hpp"
#include "Point.hpp"
#include "DiskScan.hpp"

namespace pyscan {

    bool colinear(pt2_t const& pt1,
                  pt2_t const& pt2,
                  pt2_t const& pt3){
        double x1 = pt1(0), x2 = pt2(0), x3 = pt3(0),
        y1 = pt1(1), y2 = pt2(1), y3 = pt3(1);

        double
                a11 = x2 - x1, a12 = y2 - y1,
                a21 = x2 - x3, a22 = y2 - y3;
        return (a11 * a22 - a12 * a21 == 0);
    }



    void solveCircle3(pt2_t const& pt1, pt2_t const& pt2, pt2_t const& pt3,
                      double &a, double &b) {
        double x1 = pt1(0), x2 = pt2(0), x3 = pt3(0),
            y1 = pt1(1), y2 = pt2(1), y3 = pt3(1);

        // Setup a matrix equation of the form Ax = b

        // A
        double a11 = x2 - x1, a12 = y2 - y1, a21 = x2 - x3, a22 = y2 - y3;
        // b
        double b1 = (y2 * y2 + x2 * x2 - y1 * y1 - x1 * x1) / 2, b2 = (y2 * y2
                                                                       + x2 * x2 - y3 * y3 - x3 * x3) / 2;

        double detA = a11 * a22 - a12 * a21;
        // inverse of A
        double ai11 = a22 / detA, ai12 = -a12 / detA, ai21 = -a21 / detA, ai22 = a11
                                                                                 / detA;

        // A^{-1} b = x
        a = ai11 * b1 + ai12 * b2;
        b = ai21 * b1 + ai22 * b2;
    }

    void solveCircle3(pt2_t const& pt1, pt2_t const& pt2, pt2_t const& pt3,
                     double &a, double &b, double &r) {
        double x1 = pt1(0), x2 = pt2(0), x3 = pt3(0),
                y1 = pt1(1), y2 = pt2(1), y3 = pt3(1);
        // Setup a matrix equation of the form Ax = b

        // A
        double a11 = x2 - x1,
                a12 = y2 - y1,
                a21 = x2 - x3,
                a22 = y2 - y3;
        // b
        double b1 = (y2 * y2 + x2 * x2 - y1 * y1 - x1 * x1) / 2,
                b2 = (y2 * y2 + x2 * x2 - y3 * y3 - x3 * x3) / 2;

        double detA = a11 * a22 - a12 * a21;
        // inverse of A
        double ai11 = a22 / detA,
                ai12 = -a12 / detA,
                ai21 = -a21 / detA,
                ai22 = a11 / detA;

        // A^{-1} b = x
        a = ai11 * b1 + ai12 * b2;
        b = ai21 * b1 + ai22 * b2;
        r = sqrt((x1 - a) * (x1 - a) + (y1 - b) * (y1 - b));
    }

    Disk::Disk(const pyscan::pt2_t &p1, const pyscan::pt2_t &p2, const pyscan::pt2_t &p3) {
        double a, b;
        solveCircle3(p1, p2, p3, a, b, this->radius);
        this->center = pt2_t(a, b, 1.0);
    }

    void findPerpVect(pt2_t const& pt1, pt2_t const& pt2, double* u, double* v) {
        double x1 = pt1(0), x2 = pt2(0),
                y1 = pt1(1), y2 = pt2(1);
        *u = y2 - y1;
        *v = x1 - x2;
    }


    double updateCounts(std::unordered_map<size_t, size_t>& curr_counts,
                        crescent_t const& adding, crescent_t const& removing) {
        /*
        * The adding and removing crescents define the
        */
        double update_diff = 0;
        for (auto& val_pair : adding) {
            auto curr_el = curr_counts.find(val_pair.label);
            //Start here.
            if (curr_el != curr_counts.end()) {
                curr_counts[curr_el->first] = curr_el->second + 1;
            } else {
                update_diff += val_pair.value;
                curr_counts[val_pair.label] = 1;
            }
        }
        for (auto& val_pair : removing) {
            auto curr_el = curr_counts.find(val_pair.label);
            //Start here.
            if (curr_el != curr_counts.end()) {
                if (curr_el->second <= 1) {
                    update_diff -= val_pair.value;
                    curr_counts.erase(curr_el);
                } else {
                    curr_counts[curr_el->first] = curr_el->second - 1;
                }
            } else {
                assert("The current set contains element that are being removed in the set");
            }
        }
        return update_diff;
    }

    std::ostream& operator<<(std::ostream& os, std::unordered_map<size_t, size_t> const& items) {
      for (auto ix : items) {
        os << ix.first << ":" << ix.second << ", ";
      }
      return os;
    }



    std::ostream& operator<<(std::ostream& os, crescent_t const& items) {
      for (auto ix : items) {
        os << ix.label  << ", ";
      }
      return os;
    }



    bool contains(double a, double b, double r, pt2_t const& pt) {
        double la = (pt(0) - a);
        double lb = (pt(1) - b);
        return la * la + lb * lb <= r * r;
    }




    std::tuple<Disk, double> disk_scan_simple(
            point_list_t const& net,
            wpoint_list_t const& sampleM,
            wpoint_list_t const& sampleB,
            discrepancy_func_t const& scan) {
        return scan_ranges3<Disk, pt2_t, wpt2_t>(net, sampleM, sampleB, scan);
    }

    std::tuple<Disk, double> disk_scan_simple_labels(
            point_list_t const& net,
            lpoint_list_t const& sampleM,
            lpoint_list_t const& sampleB,
            discrepancy_func_t const& scan) {
        return scan_ranges3<Disk, pt2_t, lpt2_t>(net, sampleM, sampleB, scan);
    }


    double get_order(pt2_t const& p1, pt2_t const& p2, pt2_t const& p) {
        //Create a vector between the two points
        double orthoX, orthoY;
        findPerpVect(p1, p2, &orthoX, &orthoY);
        double cX = (p1(0) + p2(0)) / 2.0;
        double cY = (p1(1) + p2(1)) / 2.0;
        double a, b;
        solveCircle3(p1, p2, p, a, b);
        return orthoX * (a - cX) + orthoY * (b - cY);
    }

    bool valid_pt(pt2_t const& p1, pt2_t const& p2, pt2_t const& p) {
        return !(p1.approx_eq(p) || p2.approx_eq(p) || colinear(p1, p2, p));
    }



    std::tuple<std::vector<double>, double> compute_delta(
            wpoint_list_t::const_iterator it_b,
            wpoint_list_t::const_iterator it_e,
            pt2_t const& p1,
            pt2_t const& p2,
            std::vector<double> const& orders,
            Disk const& start_disk) {

        double weight = 0;
        std::vector<double> counts(orders.size(), 0.0);
        for (; it_b != it_e; it_b++) {
            if (start_disk.contains(*it_b)) {
                weight += it_b->get_weight();
            }

            if (valid_pt(p1, p2, *it_b)) {

                auto lb = std::lower_bound(orders.begin(), orders.end(),
                                       get_order(p1, p2, *it_b));
                if (lb == orders.end()) {
                    continue;
                }
                if (start_disk.contains(*it_b)) {
                    counts[lb - orders.begin()] -= it_b->get_weight();
                } else {
                    counts[lb - orders.begin()] += it_b->get_weight();
                }
            }
        }
        return make_tuple(counts, weight);
    }

    /*
     * Computes a new set of iterators and also an ordered set of values
     */
    std::vector<double> preprocess_net(pt2_t p1, pt2_t p2, double min_dist, double max_dist, point_it_t& nB, point_it_t& nE) {
        nE = std::partition(nB, nE, [&] (pt2_t const& p) {
            if (valid_pt(p1, p2, p)) {
                double a, b, r;
                solveCircle3(p1, p2, p, a, b, r);
                return min_dist <= r && r <= max_dist;
            }
            return false;

        });
        std::sort(nB, nE, [&](pt2_t const &pt1, pt2_t const &pt2) {
            return get_order(p1, p2, pt1) < get_order(p1, p2, pt2);
        });
        std::vector<double> orderV;
        orderV.reserve(nE - nB);
        for (auto b = nB; b != nE; b++) {
            orderV.push_back(get_order(p1, p2, *b));
        }
        return orderV;
    }




    std::tuple<Disk, double> disk_scan_restricted(
            pt2_t p1,
            pt2_t p2,
            point_list_t net,
            wpoint_list_t const& sampleM,
            wpoint_list_t const& sampleB,
            double min_dist,
            double max_dist,
            double m_Total,
            double b_Total,
            const discrepancy_func_t& scan) {


        Disk currMax;
        double maxStat = 0;
        //First remove all colinear, duplicate, points.
        auto nB = net.begin(), nE = net.end();

        if (p1.approx_eq(p2)) {
            return std::make_tuple(currMax, 0);
        }
        //This modifies nB and nE
        auto orderV = preprocess_net(p1, p2, min_dist, max_dist, nB, nE);
        if (nE - nB == 0) {
            return std::make_tuple(currMax, 0);
        }
        Disk start_disk(p1, p2, *nB);
        //Compute delta
        std::vector<double> mDelta, bDelta;
        double mCount, bCount;
        std::tie(mDelta, mCount) = compute_delta(sampleM.begin(), sampleM.end(), p1, p2, orderV, start_disk);
        std::tie(bDelta, bCount) = compute_delta(sampleB.begin(), sampleB.end(), p1, p2, orderV, start_disk);

        mDelta[0] = 0, bDelta[0] = 0;

        //Now scan over the counts.
        auto size = nE - nB;
        for (int k = 0; k < size; k++) {
            mCount += mDelta[k];
            bCount += bDelta[k];
            double m_hat = mCount / m_Total;
            double b_hat = bCount / b_Total;
            double newStat = scan(m_hat, b_hat);
            if (maxStat <= newStat) {
                currMax = Disk(p1, p2, *(nB + k));
                maxStat = newStat;
            }
        }
        return std::make_tuple(currMax, maxStat);
    }


    template <typename P, typename T>
    void insert_queue(std::vector<P> &items, T b_it, T e_it) {
        for (auto b = b_it; b != e_it; b++) {
            items.push_back(b->second);
        }
    }


    std::tuple<std::vector<crescent_t>, std::vector<crescent_t>, std::unordered_map<size_t, size_t>, double>
    compute_delta(
            lpoint_list_t::const_iterator it_b,
            lpoint_list_t::const_iterator it_e,
            pt2_t const& p1,
            pt2_t const& p2,
            std::vector<double> const& orders,
            Disk const& start_disk) {

        double weight = 0;
        std::vector<crescent_t> add_counts(orders.size(), crescent_t());
        std::vector<crescent_t> remove_counts(orders.size(), crescent_t());

        std::unordered_map<size_t, size_t> labels;
        for (; it_b != it_e; it_b++) {
            if (valid_pt(p1, p2, *it_b)) {
                auto lb = std::lower_bound(orders.begin(), orders.end(),
                                           get_order(p1, p2, *it_b));
                if (lb == orders.end()) {
                    continue;
                }
                if (start_disk.contains(*it_b)) {
                    remove_counts[lb - orders.begin()].emplace_back(it_b->get_label(), it_b->get_weight());

                    if (labels.find(it_b->get_label()) == labels.end()) {
                        weight += it_b->get_weight();
                        labels.emplace(it_b->get_label(), 1);
                    } else {
                        labels[it_b->get_label()] += 1;
                    }
                } else {
                    add_counts[lb - orders.begin()].emplace_back(it_b->get_label(), it_b->get_weight());
                }
            } else {
                if (start_disk.contains(*it_b)) {

                    if (labels.find(it_b->get_label()) == labels.end()) {
                        weight += it_b->get_weight();
                        labels.emplace(it_b->get_label(), 1);
                    } else {
                        labels[it_b->get_label()] += 1;
                    }
                }
            }
        }
        return make_tuple(add_counts, remove_counts, labels, weight);
    }

    auto disk_scan_restricted(
            pt2_t p1,
            pt2_t p2,
            point_list_t net,
            lpoint_list_t const& sampleM,
            lpoint_list_t const& sampleB,
            double min_dist,
            double max_dist,
            double m_Total,
            double b_Total,
            std::function<double(double, double)> const& scan)
    -> std::tuple<Disk, double> {

        Disk currMax;
        double maxStat = 0;

        if (net.size() < 3 || p1.approx_eq(p2)) {
            return std::make_tuple(currMax, maxStat);
        }
        auto nB = net.begin(), nE = net.end();

        auto orderV = preprocess_net(p1, p2, min_dist, max_dist, nB, nE);
        if (nE - nB == 0) {
            return std::make_tuple(currMax, 0);
        }

        Disk start_disk(p1, p2, *nB);
        double m_count, b_count;
        std::vector<crescent_t> mCountsR, bCountsR, mCountsA, bCountsA;
        std::unordered_map<size_t, size_t> m_curr_set, b_curr_set;
        std::tie(mCountsA, mCountsR, m_curr_set, m_count) = compute_delta(sampleM.begin(), sampleM.end(),
                p1, p2, orderV, start_disk);
        std::tie(bCountsA, bCountsR, b_curr_set, b_count) = compute_delta(sampleB.begin(), sampleB.end(),
                p1, p2, orderV, start_disk);


        /*----------------------------------------------*/
        //Now scan over the counts.

        auto size = nE - nB;
        for (int k = 0; k < size; k++) {
            m_count += updateCounts(m_curr_set, mCountsA[k], mCountsR[k]);
            b_count += updateCounts(b_curr_set, bCountsA[k], bCountsR[k]);

            double m_hat = m_count / m_Total;
            double b_hat = b_count / b_Total;
            double newStat = scan(m_hat, b_hat);
            if (maxStat <= newStat) {
                double a, b, r;
                solveCircle3(p1, p2, *(nB + k), a, b, r);
                Disk currDisk(a, b, r);
                currMax = currDisk;
                maxStat = newStat;
            }
        }
        return std::make_tuple(currMax, maxStat);
    }

    double computeLTotal(std::vector<LPoint<>> const& pts) {
        return computeLabelTotal(pts.begin(), pts.end());
    }

    double computeLTotal(std::vector<WPoint<>> const& pts) {
        return computeTotal(pts.begin(), pts.end());
    }


    template <typename T>
    std::tuple<Disk, double> disk_scan_internal(point_list_t const& net,
                                       T const& sampleM,
                                       T const& sampleB,
                                       discrepancy_func_t const& scan) {
        double m_Total = computeLTotal(sampleM);
        double b_Total = computeLTotal(sampleB);
        Disk currMax;
        double maxStat = 0;

        for (auto pt1 = net.begin(); pt1 != net.end() - 1; pt1++) {
            for (auto pt2 = pt1 + 1; pt2 != net.end(); pt2++) {
                Disk d1;
                double possible_max;
                std::tie(d1, possible_max) = disk_scan_restricted(*pt1, *pt2,
                                                                  net,
                                                                  sampleM,
                                                                  sampleB,
                                                                  0.0,
                                                                  std::numeric_limits<double>::infinity(),
                                                                  m_Total,
                                                                  b_Total,
                                                                  scan);
                if (possible_max > maxStat) {
                    currMax = d1;
                    maxStat = possible_max;
                }
            }
        }
        return std::make_tuple(currMax, maxStat);
    }

    std::tuple<Disk, double> disk_scan(point_list_t const& net,
                                                wpoint_list_t const& sampleM,
                                                wpoint_list_t const& sampleB,
                                                discrepancy_func_t const& scan) {
        return disk_scan_internal(net, sampleM, sampleB, scan);
    }

    std::tuple<Disk, double> disk_scan_labels(point_list_t const& net,
                lpoint_list_t const& sampleM,
                lpoint_list_t const& sampleB,
                discrepancy_func_t const& scan) {
        return disk_scan_internal(net, sampleM, sampleB, scan);
    }

    template <typename T>
    std::tuple<Disk, double> disk_scan_scale(point_list_t const& net,
                                             std::vector<T> const& sampleM,
                                             std::vector<T> const& sampleB,
                                             uint32_t grid_r,
                                             discrepancy_func_t const& scan) {

        //Calculate the total measured and baseline value.

        double m_Total = computeLTotal(sampleM);
        double b_Total = computeLTotal(sampleB);


        using Pt_t = T;

        //using sp_grid = std::unordered_map<uint64_t,
        SparseGrid<Point<>> grid_net(net, grid_r);

        SparseGrid<Pt_t> grid_sample_m(sampleM, grid_r);
        SparseGrid<Pt_t> grid_sample_b(sampleB, grid_r);

        using net_it = typename SparseGrid<Point<>>::buck_it;
        using pt_it = typename SparseGrid<Pt_t>::buck_it;
        Disk currMax;
        double maxStat = 0;


        auto compareFirst = [](auto const& lhs, auto const& rhs) {
            return lhs.first < rhs.first;
        };

        for (auto center_cell = grid_net.begin(); center_cell != grid_net.end();
                    center_cell = std::upper_bound(center_cell, grid_net.end(), *center_cell, compareFirst)) {
            //Returns the next cell

            std::vector<Point<>> net_chunk;
            std::vector<Pt_t> m_sample_chunk;
            std::vector<Pt_t> b_sample_chunk;

            int32_t i, j;
            std::tie(i, j) = grid_net.get_cell(center_cell->second);
            //std::cout << i << " " << j << std::endl;
            net_it pt_begin, pt_end;
            std::tie(pt_begin, pt_end) = grid_net(i, j);

            for (int k = i - 3; k <= i + 3; k++) {// Through y, but ignore the last part
                for (int l = j - 3; l <= j + 3; l++) { // through x
                    //Rotate this queue
                    net_it chunk_begin, chunk_end;
                    std::tie(chunk_begin, chunk_end) = grid_net(k, l);
                    insert_queue(net_chunk, chunk_begin, chunk_end);

                    pt_it wchunk_begin, wchunk_end;
                    std::tie(wchunk_begin, wchunk_end) = grid_sample_m(k, l);
                    insert_queue(m_sample_chunk, wchunk_begin, wchunk_end);

                    std::tie(wchunk_begin, wchunk_end) = grid_sample_b(k, l);
                    insert_queue(b_sample_chunk, wchunk_begin, wchunk_end);
                }
            }


            if (net_chunk.size() < 3) {
                continue;
            }
            for (auto pt1 = pt_begin; pt1 != pt_end; pt1++) {
                for (auto pt2 = net_chunk.begin(); pt2 != net_chunk.end(); pt2++) {
                    if (pt1->second.approx_eq(*pt2)) {
                        continue;
                    }
                    Disk d1;
                    std::vector<Point<>> new_net = net_chunk;
                    double possible_max;
                    std::tie(d1, possible_max) = disk_scan_restricted(pt1->second, *pt2,
                                                                      new_net,
                                                                      m_sample_chunk,
                                                                      b_sample_chunk,
                                                                      grid_net.get_resolution(),
                                                                      2 * grid_net.get_resolution(),
                                                                      m_Total,
                                                                      b_Total,
                                                                      scan);
                    if (possible_max > maxStat) {
                        currMax = d1;
                        maxStat = possible_max;
                    }
                }
            }
        }
        //std::cout << maxStat << std::endl;
        return std::make_tuple(currMax, maxStat);
    }


    template<typename T>
    std::tuple<Disk, double> cached_disk_scan(
            point_list_t const& net,
            T const& sampleM,
            T const& sampleB,
            discrepancy_func_t const& scan) {
        /*
         * Computes all resolutions
         */
        Disk currMax;
        double maxStat = -std::numeric_limits<double>::infinity();
        for (uint32_t resolution = 2; resolution < 31; resolution++) {
            uint32_t grid_r = 1u << resolution;
            std::cout << grid_r << std::endl;
            auto new_max = disk_scan_scale(net, sampleM, sampleB, grid_r, scan);
            if (maxStat < std::get<1>(new_max)) {
                std::tie(currMax, maxStat) = new_max;
            }
        }
        return std::make_tuple(currMax, maxStat);
    }


    template std::tuple<Disk, double> cached_disk_scan<wpoint_list_t>(
            point_list_t const&,
            wpoint_list_t const&,
            wpoint_list_t const&,
            discrepancy_func_t const& );

    template std::tuple<Disk, double> cached_disk_scan<lpoint_list_t>(
            point_list_t const&,
            lpoint_list_t const&,
            lpoint_list_t const& ,
            discrepancy_func_t const&);
}
