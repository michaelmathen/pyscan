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
    auto randomLPointsUnique(int test_size) -> std::vector<pyscan::LPoint<>>;

    auto removeLabels(std::vector<pyscan::LPoint<>> const& pts) -> std::vector<pyscan::Point<>>;

    template <typename F>
    auto maxVecD(std::vector<pyscan::VecD> const& vecs, F func)
    -> pyscan::VecD {
      pyscan::VecD max_pt(0.0, 0.0);
      for (auto pt : vecs) {
          if (func(pt) > func(max_pt)) {
              max_pt = pt;
          }
      }
      return max_pt;
    }

}
#endif //PYSCAN_UTILITIES_HPP
