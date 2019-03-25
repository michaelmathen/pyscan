//
// Created by mmath on 3/24/19.
//

#ifndef PYSCAN_JEFFCODES_HPP
#define PYSCAN_JEFFCODES_HPP



#include "RectangleScan.hpp"

namespace pyscan {

    enum class F_Type {
        POISSON,
        GAUSSIAN,
        BERNOULLI,
        GAMMA
    };

    std::tuple<Rectangle, double> naive_find_rect(
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            F_Type func);

    std::tuple<Rectangle, double> naive_scan_grid(
            int grid_size,
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            F_Type func);

    std::tuple<Rectangle, double> scan_grid(
            int grid_size,
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            double eps,
            F_Type func);

    std::tuple<Rectangle, double> naive_approx_find_rect(
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            double eps,
            F_Type func);

    std::tuple<Rectangle, double> find_rect(
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            double eps,
            F_Type func);
};

#endif //PYSCAN_JEFFCODES_HPP
