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
}
#endif //PYSCAN_UTILITIES_HPP
