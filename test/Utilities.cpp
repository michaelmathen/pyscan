//
// Created by mmath on 10/2/17.
//

#include <random>

#include "Utilities.hpp"

namespace pyscantest {
  auto randomVec(int test_size) -> std::vector<pyscan::VecD> {
      std::random_device rd;
      std::default_random_engine generator(rd());
      std::uniform_real_distribution<double> distribution (0.0,1.0);
      std::vector<pyscan::VecD> points;
      for (int i = 0; i < test_size; i++) {
          points.emplace_back(distribution(generator), distribution(generator));
      }
      return points;
  }


    auto randomPoints(int test_size) -> std::vector<pyscan::Point<>> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);
        std::vector<pyscan::Point<>> points;
        for (int i = 0; i < test_size; i++) {
            points.emplace_back(distribution(generator), distribution(generator),
                                distribution(generator), distribution(generator));
        }
        return points;
    }

    auto randomLPoints(int test_size, size_t num_labels) -> std::vector<pyscan::LPoint<>> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);
        std::uniform_int_distribution<size_t> label_dist(0, num_labels);
        std::vector<pyscan::LPoint<>> points;
        for (int i = 0; i < test_size; i++) {
            points.emplace_back(label_dist(generator),
              distribution(generator),
              distribution(generator),
              distribution(generator),
              distribution(generator));
        }
        return points;
    }


}
