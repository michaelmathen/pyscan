//
// Created by mmath on 10/2/17.
//

#ifndef PYSCAN_TEST_UTILITIES_HPP
#define PYSCAN_TEST_UTILITIES_HPP

#include <functional>

#include "../src/Point.hpp"


namespace pyscantest {

	using Vec2 = std::array<double, 2>;

    auto randomPoints(int test_size) -> std::vector<pyscan::Point<>>;

    auto randomVec(int test_size) -> std::vector<Vec2>;

    Vec2 maxVec2(std::vector<Vec2> const& vec, 
    	std::function<double(Vec2)> const& f);

	auto randomWPoints(int test_size) -> std::vector<pyscan::WPoint<>>;
    auto randomLabels(int test_size, size_t num_labels) -> pyscan::label_list;

}
#endif //PYSCAN_UTILITIES_HPP
