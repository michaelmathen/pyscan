//
// Created by mmath on 10/2/17.
//

#include "Range.hpp"
#include "HalfSpaceScan.hpp"
#include "Test_Utilities.hpp"
#include "DiskScan.hpp"

#include <ctime>
#include <limits.h>
#include <random>
#include <iostream>

pyscan::discrepancy_func_t stat = [](double m, double m_total, double b, double b_total) {
    return fabs(m / m_total - b / b_total);
};

int main() {

    const static int n_size = 100;
    const static int s_size = 10000;
    double min_res = .1;
    double alpha = .01;
    auto n_pts = pyscantest::randomLPoints2(n_size, 20);
    auto m_pts = pyscantest::randomLPoints2(s_size, 20);
    auto b_pts = pyscantest::randomLPoints2(s_size, 20);

    auto begin = std::clock();
    auto [d1, d1value] = pyscan::max_disk_scale_labeled(n_pts, m_pts, b_pts, alpha, min_res, stat);

    auto end = std::clock();

    std::cout << d1value << " " << std::endl;
    std::cout << static_cast<double>(end - begin) / CLOCKS_PER_SEC << std::endl;

//    begin = std::clock();
//    std::tie(d1, d1value) = pyscan::max_disk_lift_labeled(n_pts, m_pts, b_pts, stat);
//
//    end = std::clock();
//
//    std::cout << d1value << " " << std::endl;
//    std::cout << static_cast<double>(end - begin) / CLOCKS_PER_SEC << std::endl;
    return 0;
}