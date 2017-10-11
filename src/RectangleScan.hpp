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

        bool contains(Point<> const& pt) {
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
        Grid(point_it begin, point_it end, size_t r_arg, bool run_net);
        Grid(point_it begin, point_it end, size_t r_arg);
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



    Rectangle maxLabeledRectStat(std::vector<LPoint<>> const& net,
                                 std::vector<LPoint<>> const& m_points,
                                 std::vector<LPoint<>> const& b_points, double rho);

    Subgrid maxSubgridKullSlow(Grid const &grid, double rho);
    Subgrid maxSubgridLinearSlow(Grid const& grid, double a, double b);

    /*
     * Simple 1/eps^3 algorithm that computes the max subgrid over a linear function.
     */
    Subgrid maxSubgridLinearSimple(Grid const &grid, double a, double b);
    Subgrid maxSubgridLinearSimpleStat(Grid const& grid, double rho);


    Subgrid maxSubgridLinKull(Grid const& grid, double eps, double rho);
    Subgrid maxSubgridLinGamma(Grid const& grid, double eps, double rho);

    Subgrid maxSubgridLinearG(Grid const &grid, size_t r_prime, double a, double b);

    //template<typename F>
    //subgrid maxSubgridConvex(Grid const& grid, int r_prime, F func) {

    //}
}
#endif //PYSCAN_RECTANGLESCAN_HPP
