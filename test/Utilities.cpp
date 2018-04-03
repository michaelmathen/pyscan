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
      points.emplace_back(0, 0);
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
        std::uniform_int_distribution<size_t> label_dist(0, num_labels - 1);
        std::vector<pyscan::LPoint<>> points;
        std::vector<double> m_map(num_labels, -1);
        std::vector<double> b_map(num_labels, -1);
        for (int i = 0; i < test_size; i++) {
            size_t label = label_dist(generator);
            if (m_map[label] < 0) {
              m_map[label] = distribution(generator);
              b_map[label] = distribution(generator);
            }
            points.emplace_back(label,
              m_map[label],
              b_map[label],
              distribution(generator),
              distribution(generator));
        }
        return points;
    }
auto randomLPointsUnique(int test_size) -> std::vector<pyscan::LPoint<>> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);
        std::vector<pyscan::LPoint<>> points;
        for (int i = 0; i < test_size; i++) {
            points.emplace_back(i,
              distribution(generator),
              distribution(generator),
              distribution(generator),
              distribution(generator));
        }
        return points;
    }

auto removeLabels(std::vector<pyscan::LPoint<>> const& pts) -> std::vector<pyscan::Point<>> {
  std::vector<pyscan::Point<>> new_pts;
  for (auto &pt : pts) {
    new_pts.push_back(pyscan::removeLabel(pt));
  }
  return new_pts;
}

}