//
// Created by mmath on 10/2/17.
//

#include <random>

#include "Utilities.hpp"

namespace pyscantest {

    auto randomPoints(int test_size) -> std::vector<pyscan::Point<>> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);
        std::vector<pyscan::Point<>> points;
        for (int i = 0; i < test_size; i++) {
            points.emplace_back(distribution(generator), distribution(generator), 1.0);
        }
        return points;
    }

    auto randomVec(int test_size) -> std::vector<std::array<double, 2>> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);
        std::vector<std::array<double, 2>> points;
        for (int i = 0; i < test_size; i++) {
            points.push_back(std::array<double, 2>{distribution(generator), distribution(generator)});
        }
        return points;
    }

    auto randomWPoints(int test_size) -> std::vector<pyscan::WPoint<>> {
        auto pts = randomVec(test_size);
        std::vector<pyscan::WPoint<>> wpoints(test_size, 0);
        std::transform(pts.begin(), pts.end(), wpoints.begin(), [](auto const& pt) {
            return pyscan::WPoint<>(1.0, pt[0], pt[1], 1.0);
        });
        return wpoints;
    }

    Vec2 maxVec2(std::vector<Vec2> const& vec,
        std::function<double(Vec2)> const& f) {

        Vec2 maxV = *vec.begin();
        for (auto v = vec.begin(); v != vec.end(); v++) {
            if (f(maxV) <= f(*v)) {
                maxV = *v;
            }
        }
        return maxV;
    }


    auto randomLabels(int test_size, size_t num_labels) -> pyscan::label_list {
        std::random_device rd;
        std::default_random_engine generator(rd());
        pyscan::label_list labels;
        using label_t = pyscan::label_list::value_type;

        std::uniform_int_distribution<label_t> label_dist(0, num_labels - 1);
        for (int i = 0; i < test_size; i++) {
            size_t label = label_dist(generator);
            labels.emplace_back(label);
        }
        return labels;
    }

// auto randomLPointsUnique(int test_size) -> std::vector<pyscan::LPoint<>> {
//         std::random_device rd;
//         std::default_random_engine generator(rd());
//         std::uniform_real_distribution<double> distribution (0.0,1.0);
//         std::vector<pyscan::LPoint<>> points;
//         for (int i = 0; i < test_size; i++) {
//             points.emplace_back(i,
//               distribution(generator),
//               distribution(generator),
//               distribution(generator),
//               distribution(generator));
//         }
//         return points;
//     }

// auto removeLabels(std::vector<pyscan::LPoint<>> const& pts) -> std::vector<pyscan::Point<>> {
//   std::vector<pyscan::Point<>> new_pts;
//   for (auto &pt : pts) {
//     new_pts.push_back(pyscan::removeLabel(pt));
//   }
//   return new_pts;
// }

}
