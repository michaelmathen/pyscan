//
// Created by mmath on 10/2/17.
//

#include "Range.hpp"
#include "HalfSpaceScan.hpp"
#include "Test_Utilities.hpp"

#include <ctime>
#include <limits.h>
#include <random>
#include <iostream>

pyscan::discrepancy_func_t stat = [](double m, double m_total, double b, double b_total) {
    return fabs(m / m_total - b / b_total);
};

int main() {

    const static int n_size = 500;
    const static int s_size = 80000;
    auto n_pts = pyscantest::randomPoints(n_size);
    auto m_pts = pyscantest::randomLPoints(s_size, s_size / 4);
    auto b_pts = pyscantest::randomLPoints(s_size , s_size / 4);

    auto begin = std::clock();
    auto [d1, d1value] = pyscan::max_halfplane_labeled(n_pts, m_pts, b_pts, stat);

    auto end = std::clock();

    std::cout << d1value << " " << d1.get_coords() << std::endl;
    std::cout << static_cast<double>(end - begin) / CLOCKS_PER_SEC << std::endl;
    return 0;
}