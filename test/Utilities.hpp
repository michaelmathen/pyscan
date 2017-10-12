//
// Created by mmath on 10/2/17.
//

#ifndef PYSCAN_TEST_UTILITIES_HPP
#define PYSCAN_TEST_UTILITIES_HPP
#include "../src/Vecky.hpp"
#include "../src/RectangleScan.hpp"
#include "../src/Point.hpp"
namespace pyscantest {
    auto randomVec(int test_size) -> std::vector<pyscan::VecD>;
    auto randomPoints(int test_size) -> std::vector<pyscan::Point<>>;
    auto randomLPoints(int test_size, size_t num_labels) -> std::vector<pyscan::LPoint<>>;
}
#endif //PYSCAN_UTILITIES_HPP
