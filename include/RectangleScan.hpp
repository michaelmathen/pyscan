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

#include "Range.hpp"
#include "Statistics.hpp"
#include "Point.hpp"
#include "Utilities.hpp"

//#define DEBUG

namespace pyscan {
    using point_it= std::vector<Point<> >::iterator;
    using subgrid = std::tuple<int, int, int, int, double>;


    class Subgrid {
        size_t u_x, u_y, l_x, l_y;
        double value;
    public:
        Subgrid(size_t ux, size_t uy, size_t lx, size_t ly, double val) :
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
            return u_x >= pt(0) && pt(0) >= l_x && u_y >= pt(1) && pt(1) >= l_y;
        }

        void setValue(double v) {
          value = v;
        }
        size_t lowX() const { return l_x; }
        size_t upX() const { return u_x; }
        size_t lowY() const { return l_y; }
        size_t upY() const { return u_y; }
        double fValue() const { return value; }
    };


    class Rectangle : public Range<2> {
        double u_x, u_y, l_x, l_y;
    public:
        Rectangle(const pt2_t& p1, const pt2_t& p2, const pt2_t& p3, const pt2_t& p4) {
            u_x = std::max({p1(0), p2(0), p3(0), p4(0)});
            l_x = std::min({p1(0), p2(0), p3(0), p4(0)});
            u_y = std::max({p1(1), p2(1), p3(1), p4(1)});
            l_y = std::min({p1(1), p2(1), p3(1), p4(1)});
        }

        Rectangle(double ux, double uy, double lx, double ly) : u_x(ux), u_y(uy), l_x(lx), l_y(ly) {}

        inline bool contains(const pt2_t& p1) const final {
            return (u_x >= p1(0) && p1(0) >= l_x && u_y >= p1(1) && p1(1) >= l_y);
        }

        inline bool intersects_segment(const pt2_t &p1, const pt2_t &p2) const final {
            //Check if any of points are inside of the rectangle.
            if (contains(p1) || contains(p2)) {
                return true;
            } else {
                //Check if any of the segments cross the line segment defined by p1<->p2.
                auto ur_pt = pt2_t(u_x, u_y, 1.0);
                auto ul_pt = pt2_t(l_x, u_y, 1.0);
                auto lr_pt = pt2_t(u_x, l_y, 1.0);
                auto ll_pt = pt2_t(l_x, l_y, 1.0);
                return crosses_segment(ur_pt, ul_pt, p1, p2) || crosses_segment(ul_pt, ll_pt, p1, p2) ||
                        crosses_segment(ur_pt, lr_pt, p1, p2) || crosses_segment(lr_pt, ll_pt, p1, p2);
            }
        }

        std::string toString() const {
            std::stringstream ss;
            ss << "Rectangle(" << u_x << ", " << u_y << ", " << l_x << ", " << l_y << ")";
            return ss.str();
        }
        double lowX() const { return l_x; }
        double upX() const { return u_x; }
        double lowY() const { return l_y; }
        double upY() const { return u_y; }
    };

    enum class I_Type {
        VALUE,
        MAX
    };

    class Grid {
        /*
         * Grid as defined in the SODA paper. Construction takes O(mlog r + r^2) time where
         * m = end - begin.
         */
        long r;
        std::vector<double> red_counts;
        std::vector<double> blue_counts;
        std::vector<double> x_coords;
        std::vector<double> y_coords;
        double total_red_weight = 0;
        double total_blue_weight = 0;
    public:
        Grid(size_t r_arg, point_list_t& red, weight_list_t& red_w, point_list_t& blue, weight_list_t& blue_w);
        Grid(point_list_t& net, point_list_t& red, weight_list_t& red_w, point_list_t& blue, weight_list_t& blue_w);
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


    std::tuple<Rectangle, double> max_rect_labeled(size_t r, lpoint_list_t const& m_points, lpoint_list_t const& b_points, const discrepancy_func_t& func);

    Subgrid maxSubgridLinearG(Grid const &grid, long r_prime, double a, double b);

    Subgrid maxSubgridLinearSimple(Grid const& grid, double eps, discrepancy_func_t const& f);
    Subgrid maxSubgridLinearSimple(Grid const &grid, double a, double b);
    Subgrid maxSubgridNonLinear(Grid const &grid, discrepancy_func_t const& func);
    Subgrid maxSubgridLinearTheory(Grid const& grid, double eps, discrepancy_func_t const& f);
}
#endif //PYSCAN_RECTANGLESCAN_HPP
