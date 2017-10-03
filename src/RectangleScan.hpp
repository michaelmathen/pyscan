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
    using point_it= std::vector<Point<int> >::iterator;
    using subgrid = std::tuple<int, int, int, int, double>;


    template<typename P_it, typename I_it, typename Comp>
    void quantiles(P_it begin, P_it end, I_it i_begin, I_it i_end, Comp comp) {
        //std::random_device rd;
        //std::mt19937 gen(rd());
        //std::uniform_real_distribution<double> dis(0, 1);
        int k = i_end - i_begin;
        double width = static_cast<double>(end - begin) / k;
        // Choose a start index randomly by either rounding down or up.
        int start_index = 0;
        /*
        if (dis(gen) > width  2) {
            start_index = (int) (width / 2);
        } else {
            start_index = (int) (width / 2 + 1);
        }
         */

        //Pick a set of evenly spread indices
        auto i_b = i_begin;
        for (int i = 0; i < k; i++) {
            //double prob = i * width + start_index;
            (*i_b) = i * width + start_index;
            i_b++;
        }
        std::sort(begin, end, comp);
        /*
        std::cout << "got here" << std::endl;
        std::deque<std::tuple<int, int>> tuple_stack;
        tuple_stack.push_back(std::make_tuple(0, (int)(i_end - i_begin) - 1));
        while (!tuple_stack.empty()) {
            auto &part = tuple_stack.back();
            tuple_stack.pop_back();
            int lo = std::get<0>(part);
            int hi = std::get<1>(part);
            if (lo < hi) {
                int mid = (lo + hi) / 2;
                std::cout << hi << " " << mid << " " << lo << std::endl;
                int i_lo = *(i_begin + lo);
                int i_hi = *(i_begin + hi);
                int i_mid = *(i_begin + mid);
                std::nth_element(begin + i_lo, begin + i_hi, begin + i_mid, comp);
                tuple_stack.push_back(std::make_tuple(lo, mid - 1));
                tuple_stack.push_back(std::make_tuple(mid + 1, hi));
            }
        }
         */
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

        bool contains(Point<int> const& pt) {
            return u_x >= get<0>(pt) && get<0>(pt) >= l_x && u_y >= get<1>(pt) && get<1>(pt) >= l_y;
        }

        Bound_t lowX() const { return l_x; }
        Bound_t upX() const { return u_x; }
        Bound_t lowY() const { return l_y; }
        Bound_t upY() const { return u_y; }
        double fValue() const { return value; }
    };

    using Rectangle = RectBase<double>;
    using Subgrid = RectBase<size_t>;

    template<typename Weight=int>
    class Grid {
        /*
         * Grid as defined in the SODA paper. Construction takes O(mlog r + r^2) time where
         * m = end - begin.
         */
        size_t r;
        std::vector<Weight> red_counts;
        std::vector<Weight> blue_counts;
        std::vector<double> x_coords;
        std::vector<double> y_coords;
        Weight total_red_weight = 0;
        Weight total_blue_weight = 0;
    public:

        Grid(point_it begin, point_it end, size_t r_arg, bool run_net) : r(r_arg),
                                                                  red_counts(r_arg * r_arg, 0),
                                                                  blue_counts(r_arg * r_arg, 0),
                                                                  x_coords(),
                                                                  y_coords()
        {
            std::cout << "NetRect constructor called" << std::endl;
            std::random_device rd;
            std::default_random_engine rand_eng(rd());

            //std::cout << "before the shuffle" << std::endl;
            std::shuffle (begin, end, rand_eng);
            //std::cout << "got here" << std::endl;

            size_t red_found_count = 0;
            size_t blue_found_count = 0;
            auto b = begin;

            while (red_found_count < r/2 || blue_found_count < r/2) {
                if (b->getBlueWeight() >= 1 && blue_found_count < r/2) {
                    x_coords.push_back(get<0>(*b));
                    y_coords.push_back(get<1>(*b));
                    blue_found_count += 1;
                }

                if (b->getRedWeight() >= 1 && red_found_count < r/2) {
                    x_coords.push_back(get<0>(*b));
                    y_coords.push_back(get<1>(*b));
                    red_found_count += 1;
                }
                b++;
                //std::cout << red_found_count << " " << blue_found_count << " " << r/2 << std::endl;
            }
            //std::cout << "stuck in the loop forever" << std::endl;
            std::sort(x_coords.begin(), x_coords.end());
            std::sort(y_coords.begin(), y_coords.end());
            for (auto point_it = begin; point_it != end; point_it++) {
                size_t ix = std::upper_bound(x_coords.begin(), x_coords.end(), get<0>(*point_it)) - x_coords.begin() - 1;
                size_t iy = std::upper_bound(y_coords.begin(), y_coords.end(), get<1>(*point_it)) - y_coords.begin() - 1;
                if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                    red_counts[iy * r + ix] += point_it->getRedWeight();
                    blue_counts[iy * r + ix] += point_it->getBlueWeight();
                }
                total_red_weight += point_it->getRedWeight();
                total_blue_weight += point_it->getBlueWeight();
            }
        }

        Grid(point_it begin, point_it end, size_t r_arg) : r(r_arg),
                                                           red_counts(r_arg * r_arg, 0),
                                                           blue_counts(r_arg * r_arg, 0),
                                                           x_coords(),
                                                           y_coords()
        {
            std::cout << "StandardRect constructor called" << std::endl;
            auto compX = [](Point<int> const &p1, Point<int> const &p2) {
                return get<0>(p1) < get<0>(p2);
            };
            auto compY = [](Point<int> const &p1, Point<int> const &p2) {
                return get<1>(p1) < get<1>(p2);
            };
            std::vector<int> indices(r, 0);
            quantiles(begin, end, indices.begin(), indices.end(), compX);
            x_coords.reserve(r);
            y_coords.reserve(r);
            for (auto i : indices) {
                x_coords.push_back(get<0>(*(begin + i)));
            }
            quantiles(begin, end, indices.begin(), indices.end(), compY);
            for (auto i : indices) {
                y_coords.push_back(get<1>(*(begin + i)));
            }
            std::sort(x_coords.begin(), x_coords.end());
            std::sort(y_coords.begin(), y_coords.end());
            for (auto point_it = begin; point_it != end; point_it++) {
                size_t ix = std::upper_bound(x_coords.begin(), x_coords.end(), get<0>(*point_it)) - x_coords.begin() - 1;
                size_t iy = std::upper_bound(y_coords.begin(), y_coords.end(), get<1>(*point_it)) - y_coords.begin() - 1;
                if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                    red_counts[iy * r + ix] += point_it->getRedWeight();
                    blue_counts[iy * r + ix] += point_it->getBlueWeight();
                }
                total_red_weight += point_it->getRedWeight();
                total_blue_weight += point_it->getBlueWeight();
            }
        }

        Weight totalRedWeight() const {
            return total_red_weight;
        }

        Weight totalBlueWeight() const {
            return total_blue_weight;
        }

        Weight redCount(int row, int col) const {
            return red_counts[row * r + col];
        }

        Weight blueCount(int row, int col) const {
            return blue_counts[row * r + col];
        }

        double redWeight(int row, int col) const {
            return red_counts[row * r + col] / (double)total_red_weight;
        }

        double blueWeight(int row, int col) const {
            return blue_counts[row * r + col] / (double)total_blue_weight;
        }

        double yCoord(int row) const {
            return y_coords[row];
        }

        double xCoord(int col) const {
            return x_coords[col];
        }

        size_t size() const {
            return r;
        }

        Rectangle toRectangle(Subgrid const &sg) const {
            return Rectangle(xCoord(sg.upX()), yCoord(sg.upY()),
                             xCoord(sg.lowX()), yCoord(sg.lowY()),
                             sg.fValue());
        }
    };


    Rectangle maxLabeledRectStat(std::vector<LPoint<int>> const& net,
                             std::vector<LPoint<int>> const& m_points,
                             std::vector<LPoint<int>> const& b_points);

    Rectangle maxLabeledRectStat(std::vector<LPoint<double>> const& net,
                                 std::vector<LPoint<double>> const& m_points,
                                 std::vector<LPoint<double>> const& b_points);

    Subgrid maxSubgridKullSlow(Grid<double> const &grid, double rho);
    Subgrid maxSubgridLinearSlow(Grid<int> const& grid, double a, double b);
    Subgrid maxSubgridKullSlow(Grid<int> const &grid, double rho);
    Subgrid maxSubgridLinearSlow(Grid<double> const& grid, double a, double b);

    /*
     * Simple 1/eps^3 algorithm that computes the max subgrid over a linear function.
     */
    Subgrid maxSubgridLinearSimple(Grid<double> const &grid, double a, double b);
    Subgrid maxSubgridLinearSimple(Grid<int> const &grid, double a, double b);

    Subgrid maxSubgridLinearSimpleStat(Grid<int> const& grid, double rho);

    class ValueInterval {
        size_t left;
        double value;
        size_t right;
    public:
        ValueInterval() : left(0), value(0), right(0) {}
        ValueInterval(double val, size_t left, size_t right) : left(left),
                                                               value(val),
                                                               right(right)
        {}
        void print(std::ostream& os) const {
            os << "[" << left << ", " << right << "]";
            os << " = " << value;
        }
        size_t getLeft() const { return left; }
        size_t getRight() const { return right; }
        void setLeft(size_t left_c) { left = left_c; }
        void setRight(size_t right_c) { right = right_c; }
        double getValue() const {
            return value;
        }

        void setValue(double val) {
            value = val;
        }

        ValueInterval &operator+=(double val) {
            value += val;
            return *this;
        }

    };

    class MaxInterval {
        ValueInterval left_max;
        ValueInterval right_max;
        ValueInterval center_max;
    public:
        MaxInterval() : left_max(), right_max(), center_max() {}
        MaxInterval(double value, size_t index) : left_max(value, index, index),
                                                  right_max(value, index, index),
                                                  center_max(value, index, index) {}


        MaxInterval &operator+=(MaxInterval const &op) {

            ValueInterval merged_max;
            //Figure out the new center max value
            if (right_max.getValue() < 0 || op.left_max.getValue() < 0) {
                if (right_max.getValue() < op.left_max.getValue()) {
                    merged_max = op.left_max;
                } else {
                    merged_max = right_max;
                }
            } else {
                merged_max = ValueInterval(right_max.getValue() + op.left_max.getValue(),
                                           right_max.getLeft(), op.left_max.getRight());
            }

            if (merged_max.getValue() >= center_max.getValue() && merged_max.getValue() >= op.center_max.getValue()) {
                center_max = merged_max;
            } else if (op.center_max.getValue() > center_max.getValue()) {
                center_max = op.center_max;
            } // else it remains the same.
            //Figure out the new left max value
            if (left_max.getRight() == this->getRight() && op.left_max.getValue() >= 0) {
                left_max.setRight(op.left_max.getRight());
                left_max.setValue(op.left_max.getValue() + left_max.getValue());
            }
            //Figure out the new right max value
            if (op.right_max.getLeft() == op.getLeft() && right_max.getValue() >= 0) {
                //Left boundary of the right max doesn't change
                right_max.setValue(op.right_max.getValue() + right_max.getValue());
            } else {
                // New left boundary of the right max.
                right_max = op.right_max;
            }
            right_max.setRight(op.getRight());
            return *this;
        }

        size_t getLeft() const { return left_max.getLeft(); }
        size_t getRight() const { return right_max.getRight(); }
        size_t lowCIx() const {
            return center_max.getLeft();
        }

        size_t upCIx() const {
            return center_max.getRight();
        }

        double lValue() const {
            return left_max.getValue();
        }

        double rValue() const {
            return right_max.getValue();
        }

        double cValue() const {
            return center_max.getValue();
        }

        void print(std::ostream& os) const {
            os << "[" << getLeft() << ", " << getRight() << "]";
            os << "= {";
            os << left_max << ", " << center_max << ", " << right_max << "}";
        }
    };

    enum class I_Type {
        VALUE,
        MAX
    };

    /*
    auto operator<<(std::ostream& os, const I_Type & t) -> std::ostream& {
        if (t == I_Type::VALUE) {
            os << "VALUE";
        } else {
            os << "MAX";
        }
        return os;
    }
     */

    class MaximumIntervals {
        std::vector<I_Type> intervals;
        std::vector<MaxInterval> max_intervals;
        std::vector<double> weights;
        // one to one mapping with weights
        std::vector<size_t> weight_index;

        size_t left_col;
        size_t right_col;

        MaximumIntervals(size_t l, size_t r) : left_col(l), right_col(r) {};
    public:

        void print(std::ostream& os) const {
            os << "[" << left_col << ", " << right_col << "]";
            os << max_intervals << std::endl;
            os << weights << std::endl;
            os << weight_index << std::endl;

            //os << intervals << std::endl;
        }

        size_t getWeightNum() const {
            return weights.size();
        }

        size_t getIntervalNum() const {
            return intervals.size();
        }

        MaximumIntervals(size_t r) : intervals(r, I_Type::VALUE),
                                     max_intervals(),
                                     weights(r, 0),
                                     weight_index(r, 0),
                                     left_col(0),
                                     right_col(r - 1) {
            for (size_t i = 0; i < r; i++)
                weight_index[i] = i;
        }

        void setBounds(size_t l_c, size_t r_c) {
            left_col = l_c;
            right_col = r_c;
        }

        void updateBounds(size_t l_c, size_t r_c) {
            left_col = std::min(l_c, left_col);
            right_col = std::max(r_c, right_col);
        }

        /*
         * Indices and weights are assumed to be presorted.
         */
        void updateWeights(std::vector<double> const &new_weights,
                           std::vector<size_t> const &indices,
                           double w) {

            size_t j = 0;
            for (size_t i = 0; i < new_weights.size(); i++) {
                while (j < weight_index.size() && weight_index[j] <= indices[i]) {
                    if (weight_index[j] == indices[i]) {
                        weights[j] += new_weights[i] * w;
                    }
                    j++;
                }
            }
        }

        MaximumIntervals mergeZeros(BloomFilter const& f1) const {
            auto mightBePresent = [&f1](uint32_t index) {
                return f1.mightBePresent(index);
            };

            return mergeZeros(mightBePresent);
        }
        MaximumIntervals mergeZeros(BloomFilter const& f1, BloomFilter const& f2) const {
            auto mightBePresent = [&f1, &f2](uint32_t index) {
                return f1.mightBePresent(index) || f2.mightBePresent(index);
            };
            return mergeZeros(mightBePresent);
        }

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
                            new_max.max_intervals.push_back(MaxInterval(weights[w_ix], weight_index[w_ix]));
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

        Subgrid getMax() const {
            if (intervals.size() > 0) {
                int mx_ix = 0;
                int w_ix = 0;
                MaxInterval new_max;
                if (intervals[0] == I_Type::VALUE) {
                    new_max = MaxInterval(weights[0], weight_index[0]);
                    w_ix += 1;
                } else {
                    new_max = max_intervals[0];
                    mx_ix += 1;
                }
                for (size_t i = 1; i < intervals.size(); i++) {
                    if (intervals[i] == I_Type::VALUE) {
                        new_max += MaxInterval(weights[w_ix], weight_index[w_ix]); // Merge the previous
                        w_ix += 1;
                    } else {
                        new_max += max_intervals[mx_ix]; // Merge the previous
                        mx_ix += 1;
                    }
                }
                return Subgrid(new_max.upCIx(), this->getRight(), new_max.lowCIx(), this->getLeft(), new_max.cValue());
            } else {
                return Subgrid(this->getLeft(), 0, this->getRight(), 0, -std::numeric_limits<double>::infinity());
            }
        }

        size_t maxIntervalNum() const {
            return max_intervals.size();
        }

        size_t valueNum() const {
            return weights.size();
        }

        MaxInterval getFirstMax() const {
            return max_intervals[0];
        }

        size_t getLeft() const {
            return left_col;
        }

        size_t getRight() const {
            return right_col;
        }

        void setLeft(size_t left) {
            left_col = left;
        }

        void setRight(size_t right) {
            right_col = right;
        }
    };


    struct SlabNode {
        std::vector<double> red_weights;
        std::vector<double> blue_weights;
        std::vector<size_t> indices;
        BloomFilter non_zeros;

        size_t left_col;
        size_t right_col;
        using slab_ptr = std::shared_ptr<SlabNode>;
        slab_ptr left = nullptr;
        slab_ptr right = nullptr;


        template<typename Weight>
        SlabNode(Grid<Weight> const& g, size_t left_col, size_t right_col, size_t r_prime) :
                left_col(left_col),
                right_col(right_col) {
            double red_weight = 0;
            double blue_weight = 0;
            //std::cout << "[" << left_col << ", " << right_col << "]" << std::endl;
            for (size_t row = 0; row < g.size(); row++) {
                for (size_t col = left_col; col <= right_col; col++) {
                    red_weight += g.redWeight(row, col);
                    blue_weight += g.blueWeight(row, col);
                }
                if (red_weight + blue_weight > 1 / static_cast<double>(r_prime)) {
                    red_weights.push_back(red_weight);
                    blue_weights.push_back(blue_weight);
                    red_weight = 0;
                    blue_weight = 0;
                    indices.push_back(row);
                }
            }
            //std::cout << red_weights << std::endl;
            //std::cout << indices << std::endl;
        }

        size_t getMid() const {
            return (left_col + right_col) / 2;
        }

        /*
        double measureInterval(size_t l_r, size_t u_r, double a, double b) const {
           double sum = 0;
           //std::cout << l_r << " " << u_r << std::endl;
           //std::cout << indices << std::endl;
           for (size_t i = 0; i < indices.size(); i++) {
               if (l_r <= indices[i] && indices[i] <= u_r) {
                   sum += a * red_weights[i] + b * blue_weights[i];
               }
            }
            return sum;
        }


        double measureRect(size_t l_c, size_t l_r, size_t u_c, size_t u_r, double a, double b) const {
            // If it completely overlaps
            auto measureRectH = [&](slab_ptr ptr) {
                if (ptr != nullptr) {
                    return ptr->measureRect(l_c, l_r, u_c, u_r, a, b);
                } else {
                    return 0.0;
                }
            };
            //std::cout << l_c << " " << l_r << " " << u_c << " " << u_r << std::endl;
            //std::cout << left_col << " " << right_col << std::endl;
            //print(std::cout);
            //std::cout << std::endl;

            if (l_c <= left_col && right_col <= u_c) {
                return measureInterval(l_r, u_r, a, b);
            } else if (right_col < l_c || u_c < left_col) {
                return 0.0;
            } else {
                return measureRectH(left) + measureRectH(right);
            }
        }
        */

        bool hasChildren() const {
            return left != nullptr;
        }

        void print(std::ostream& os) const {
            os << "red= " << red_weights << std::endl;
            os << "blue= " << blue_weights << std::endl;
            os << "indices= " << indices << std::endl;
        }

    };


    using slab_ptr = SlabNode::slab_ptr;

    template<typename Weight=bool>
    struct SlabTree {
        size_t r;
        slab_ptr root;
        /*
         * red_b, red_e -- iterators to the red points
         * blue_b, blue_e --iterators to the blue points
         * r -- the number of points to use in the top level slab.
         */
        SlabTree(Grid<Weight> const &grid, size_t r_prime) : r(grid.size()), root(nullptr) {
            size_t upper = r - 1, lower = 0;
            std::deque<std::tuple<int, int, slab_ptr>> stack_nodes;
            root = std::make_shared<SlabNode>(grid, lower, upper, r_prime);
            stack_nodes.push_back(std::make_tuple(lower, upper, root));
            while (stack_nodes.size() > 0) {
                auto node = stack_nodes.back();
                if (std::get<1>(node) <= std::get<0>(node)) {
                    stack_nodes.pop_back();
                    continue;
                }

                auto mid = (std::get<1>(node) + std::get<0>(node)) / 2;

                auto &parent = std::get<2>(node);
                auto left_child = std::make_shared<SlabNode>(grid, std::get<0>(node), mid, r_prime);
                auto right_child = std::make_shared<SlabNode>(grid, mid + 1, std::get<1>(node), r_prime);
                parent->left = left_child;
                parent->right = right_child;
                stack_nodes.pop_back();
                stack_nodes.push_back(std::make_tuple(std::get<0>(node), mid, left_child));
                stack_nodes.push_back(std::make_tuple(mid + 1, std::get<1>(node), right_child));
            }

            computeBlooms(root);
        }

        void print(std::ostream& os) const {
            printHelper(root, os);
        }
        /*
        double measure(size_t l_c, size_t l_r, size_t u_c, size_t u_r, double a, double b) const {
            if (root == nullptr) {
                return 0.0;
            } else {
                return root->measureRect(l_c, l_r, u_c, u_r, a, b);
            }
        }

        double measureSubgrid(Subgrid const& subgrid, double a, double b) const {
            return measure(subgrid.upX(), subgrid.lowRow(), subgrid.upCol(), subgrid.upRow(), a, b);
        }
         */

    private:

        std::vector<int> computeBlooms(slab_ptr &curr_root) {
#ifdef DEBUG
            std::cout << "Computing blooms in slab tree" << std::endl;
#endif
            if (curr_root == nullptr) {
                return std::vector<int>();
            } else {
                //curr_root->print(std::cout);
                std::vector<int> tmp_indices;
                std::vector<int> new_indices;
                auto left_non_zeros = computeBlooms(curr_root->left);
                auto right_non_zeros = computeBlooms(curr_root->right);
                //std::cout << "left = " << left_non_zeros << std::endl;
                // std::cout << "right = " << right_non_zeros << std::endl;
                auto &curr_indices = curr_root->indices;
                //std::cout << "current = " << curr_root->indices << std::endl;
                std::set_union(left_non_zeros.begin(), left_non_zeros.end(),
                               right_non_zeros.begin(), right_non_zeros.end(),
                               std::back_inserter(tmp_indices)
                );
                std::set_union(tmp_indices.begin(), tmp_indices.end(),
                               curr_indices.begin(), curr_indices.end(),
                               std::back_inserter(new_indices));
                //std::cout << "merged = " << new_indices << std::endl;
                curr_root->non_zeros = BloomFilter(new_indices.size(), .1);

                for (size_t i = 0; i < new_indices.size(); i++) {
                    //std::cout << curr_root->non_zeros.mightBePresent(i) << std::endl;
                    curr_root->non_zeros.insert(new_indices[i]);
                }
                return new_indices;
            }
        }

        void printHelper(slab_ptr curr_root, std::ostream& os) const {
            if (curr_root != nullptr) {
                curr_root->print(os);
                printHelper(curr_root->left, os);
                printHelper(curr_root->right, os);
            }
        }
    };


    Subgrid maxSubgridLinearOneSided(slab_ptr left, slab_ptr right, bool left_side, MaximumIntervals const& maxInt, double a, double b);
    Subgrid maxSubgridSpan(slab_ptr left, slab_ptr right, MaximumIntervals const& maxInt,
                         double a, double b);
    Subgrid maxSubgridLinear(slab_ptr slab, MaximumIntervals const& maxInt, double a, double b);

    template<typename Weight>
    Subgrid maxSubgridLinearG(Grid<Weight> const &grid, size_t r_prime, double a, double b) {
        SlabTree<Weight> slabT(grid, r_prime);
#ifdef DEBUG
        std::cout<<"finished the slab tree" << std::endl;
#endif
        //std::cout << slabT << std::endl;
        MaximumIntervals maxInt(grid.size());
        return maxSubgridLinear(slabT.root, maxInt, a, b);
    }

    //template<typename F>
    //subgrid maxSubgridConvex(Grid const& grid, int r_prime, F func) {

    //}
}
#endif //PYSCAN_RECTANGLESCAN_HPP
