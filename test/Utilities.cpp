//
// Created by mmath on 10/2/17.
//

#include <random>

#include "Utilities.hpp"

namespace pyscantest {
  auto random_Vec(int test_size) -> std::vector<pyscan::VecD> {
      std::random_device rd;
      std::default_random_engine generator(rd());
      std::uniform_real_distribution<double> distribution (0.0,1.0);
      std::vector<pyscan::VecD> points;
      for (int i = 0; i < test_size; i++) {
          points.emplace_back(distribution(generator), distribution(generator));
      }
      return points;
  }
}
