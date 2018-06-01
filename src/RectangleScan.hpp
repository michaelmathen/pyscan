//
// Created by mmath on 5/28/17.
//

#ifndef PYSCAN_RECTANGLESCAN_HPP
#define PYSCAN_RECTANGLESCAN_HPP
#include <memory>
#include <vector>
#include <forward_list>
#include <algorithm>
#include <cmath>
#include <random>
#include <deque>
#include <functional>
#include <tuple>
#include <limits>
#include <iostream>
#include "BloomFilter.hpp"


#include "Statistics.hpp"
#include "Point.hpp"
#include "Utilities.hpp"

//#define DEBUG

namespace pyscan {
    using point_it= std::vector<Point<> >::iterator;
    using subgrid = std::tuple<int, int, int, int, double>;


    template<typename P_it, typename I_it, typename Comp, typename Wf>
    void quantiles(P_it begin, P_it end, I_it i_begin, I_it i_end, Comp comp, Wf wf) {
        std::sort(begin, end, comp);

        //std::random_device rd;
        //std::mt19937 gen(rd());
        //std::uniform_real_distribution<double> dis(0, 1);
        int r = i_end - i_begin;
        double total_weight = 0;
        std::for_each(begin, end, [&](decltype(*begin) const& pt) {
            total_weight += wf(pt);
        });

        double eps_s = total_weight / r;
        double curr_weight = 0;
        std::for_each(begin, end, [&](decltype(*begin) const& pt) {
            curr_weight += wf(pt);
            if (curr_weight > eps_s) {
                *i_begin = pt;
                i_begin++;
                curr_weight = 0;
            }
        });
    }

    template<typename Bound_t>
    class RectBase {
        Bound_t u_x, u_y, l_x, l_y;
        double value;
    public:
        RectBase(Bound_t ux, Bound_t uy, Bound_t lx, Bound_t ly, double val) :
                u_x(ux),
                u_y(uy),
                l_x(lx),
                l_y(ly),
                value(val) {}


        std::string toString() const {
           std::stringstream ss;
           ss << "(" << u_x << " " << u_y << " " << l_x << " " << l_y << " " << value << ")";
           return ss.str();
        }

        void print(std::ostream& os) const {
           os <<  "(" << u_x << " " << u_y << " " << l_x << " " << l_y << " " << value << ")";
        }

        bool contains(Point<> const& pt) const {
            return u_x >= get<0>(pt) && get<0>(pt) >= l_x && u_y >= get<1>(pt) && get<1>(pt) >= l_y;
        }

        void setValue(double v) {
          value = v;
        }
        Bound_t lowX() const { return l_x; }
        Bound_t upX() const { return u_x; }
        Bound_t lowY() const { return l_y; }
        Bound_t upY() const { return u_y; }
        double fValue() const { return value; }
    };

    using Rectangle = RectBase<double>;
    using Subgrid = RectBase<size_t>;

    enum class I_Type {
        VALUE,
        MAX
    };

    class Grid {
        /*
         * Grid as defined in the SODA paper. Construction takes O(mlog r + r^2) time where
         * m = end - begin.
         */
        size_t r;
        std::vector<double> red_counts;
        std::vector<double> blue_counts;
        std::vector<double> x_coords;
        std::vector<double> y_coords;
        double total_red_weight = 0;
        double total_blue_weight = 0;
    public:
        Grid(size_t r_arg, point_it r_begin, point_it r_end, point_it b_begin, point_it b_end);
        Grid(point_it net_begin, point_it net_end, point_it r_begin, point_it r_end, point_it b_begin, point_it b_end);
        double totalRedWeight() const;
        double totalBlueWeight() const;
        double redCount(size_t row, size_t col) const;
        double blueCount(size_t row, size_t col) const;
        double redWeight(size_t row, size_t col) const;
        double blueWeight(size_t row, size_t col) const;
        double redSubWeight(Subgrid const& sg) const;
        double blueSubWeight(Subgrid const& sg) const;

        double yCoord(size_t row) const;
        double xCoord(size_t col) const;

        size_t size() const;
        Rectangle toRectangle(Subgrid const &sg) const;
    };

    class ValueInterval {
        size_t left;
        double value;
        size_t right;
    public:
        ValueInterval();
        ValueInterval(double val, size_t left, size_t right);
        void print(std::ostream& os) const;
        size_t getLeft() const;
        size_t getRight() const;
        void setLeft(size_t left_c);
        void setRight(size_t right_c);
        double getValue() const;
        void setValue(double val);
        ValueInterval &operator+=(double val);
    };


    class MaxInterval {
        ValueInterval left_max;
        ValueInterval right_max;
        ValueInterval center_max;
    public:
        MaxInterval();
        MaxInterval(double value, size_t index);
        MaxInterval &operator+=(MaxInterval const &op);
        size_t getLeft() const;
        size_t getRight() const;
        size_t lowCIx() const;
        size_t upCIx() const;
        double lValue() const;
        double rValue() const;
        double cValue() const;
        void print(std::ostream& os) const;
    };

    struct SlabNode {
        std::vector<double> red_weights;
        std::vector<double> blue_weights;
        std::vector<size_t> indices;
        BloomFilter non_zeros;

        using slab_ptr = std::shared_ptr<SlabNode>;
        size_t left_col;
        size_t right_col;
        slab_ptr left = nullptr;
        slab_ptr right = nullptr;


        SlabNode(Grid const &g, size_t left_col, size_t right_col, size_t r_prime);

        size_t getMid() const;

        bool hasChildren() const;

        void print(std::ostream &os) const;

    };

    class MaximumIntervals {
        std::vector<I_Type> intervals;
        std::vector<MaxInterval> max_intervals;
        std::vector<double> weights;
        // one to one mapping with weights
        std::vector<size_t> weight_index;

        size_t left_col;
        size_t right_col;

        MaximumIntervals(size_t l, size_t r);
    public:

        void print(std::ostream& os) const;

        size_t getWeightNum() const;
        size_t getIntervalNum() const;
        explicit MaximumIntervals(size_t r);
        void setBounds(size_t l_c, size_t r_c);
        void updateBounds(size_t l_c, size_t r_c);

        /*
         * Indices and weights are assumed to be presorted.
         */
        void updateWeights(std::vector<double> const &new_weights,
                           std::vector<size_t> const &indices,
                           double w);

        MaximumIntervals mergeZeros(BloomFilter const& f1) const;
        MaximumIntervals mergeZeros(BloomFilter const& f1, BloomFilter const& f2) const;
        template<typename F>
        MaximumIntervals mergeZeros(F const& mightBePresent) const {
            MaximumIntervals new_max(left_col, right_col);
            new_max.max_intervals.reserve(max_intervals.size() / 2);
            new_max.weights.reserve(weights.size() / 2);
            new_max.intervals.reserve(max_intervals.size() / 2);
            int mx_ix = 0;
            int w_ix = 0;
            /*
             std::vector<int> zeros;
            for (int i = 0; i < 20; i++) {
                if (!mightBePresent(i)) {
                    zeros.push_back(i);
                }
            }
            */
            //std::cout << zeros << std::endl;
            //std::cout << *this << std::endl;

            //std::cout << bloom << std::endl;
            bool prev_Max = false;
            for (size_t i = 0; i < intervals.size(); i++) {
                if (intervals[i] == I_Type::VALUE) {
                    if (!mightBePresent((uint32_t) weight_index[w_ix])) {
                        if (prev_Max) {
                            auto prev_el = new_max.max_intervals.end() - 1;
                            (*prev_el) += MaxInterval(weights[w_ix], weight_index[w_ix]); // Merge the previous
                        } else {
                            new_max.max_intervals.emplace_back(weights[w_ix], weight_index[w_ix]);
                            new_max.intervals.push_back(I_Type::MAX);
                        }
                        prev_Max = true;
                    } else {
                        prev_Max = false;
                        new_max.weights.push_back(weights[w_ix]);
                        new_max.weight_index.push_back(weight_index[w_ix]);
                        new_max.intervals.push_back(I_Type::VALUE);
                    }
                    w_ix += 1;
                } else {
                    if (prev_Max) {
                        auto prev_el = new_max.max_intervals.end() - 1;
                        (*prev_el) += max_intervals[mx_ix]; // Merge the previous
                    } else {
                        new_max.max_intervals.push_back(max_intervals[mx_ix]);
                        new_max.intervals.push_back(I_Type::MAX);
                    }
                    prev_Max = true;
                    mx_ix += 1;
                }
            }
            //std::cout << new_max << std::endl;
            //std::cout << std::endl;
            return new_max;
        }

        Subgrid getMax() const;
        size_t maxIntervalNum() const;
        size_t valueNum() const;
        MaxInterval getFirstMax() const;
        size_t getLeft() const;
        size_t getRight() const;
        void setLeft(size_t left);
        void setRight(size_t right);
    };

    Rectangle maxRectStatLabels(std::vector<LPoint<>> const& net,
                                 std::vector<LPoint<>> const& m_points,
                                 std::vector<LPoint<>> const& b_points,
                                 double rho);

    Subgrid maxSubgridKullSlow(Grid const &grid, double rho);
    Subgrid maxSubgridLinearSlow(Grid const& grid, double a, double b);

    /*
     * Simple 1/eps^3 algorithm that computes the max subgrid over a linear function.
     */
    Subgrid maxSubgridLinearSimple(Grid const &grid, double a, double b);
    Subgrid maxSubgridLinearSimpleStat(Grid const& grid);


    Subgrid maxSubgridLinKull(Grid const& grid, double eps, double rho);
    Subgrid maxSubgridLinKullTheory(Grid const& grid, double eps, double rho);
    Subgrid maxSubgridLinGamma(Grid const& grid, double eps);

    Subgrid maxSubgridLinearG(Grid const &grid, int r_prime, double a, double b);

    Rectangle maxRectSlowStatLabels(std::vector<LPoint<>> const& net,
                                 std::vector<LPoint<>> const& m_points,
                                 std::vector<LPoint<>> const& b_points,
                                 double rho);
    //template<typename F>
    //subgrid maxSubgridConvex(Grid const& grid, int r_prime, F func) {

    //}
}
#endif //PYSCAN_RECTANGLESCAN_HPP
