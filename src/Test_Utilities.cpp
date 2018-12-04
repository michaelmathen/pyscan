//
// Created by mmath on 10/2/17.
//

#include <random>
#include <algorithm>
#include <Point.hpp>


#include "Test_Utilities.hpp"

namespace pyscantest {

    auto randomPoints2(int test_size) -> std::vector<pyscan::pt2_t> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);
        std::vector<pyscan::pt2_t> points;
        for (int i = 0; i < test_size; i++) {
            points.emplace_back(distribution(generator), distribution(generator), 1.0);
        }
        return points;
    }

    auto randomLabels(size_t test_size, size_t num_labels) -> pyscan::label_list_t {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::vector<size_t> labels;
        using label_t = size_t;

        std::uniform_int_distribution<label_t> label_dist(0, num_labels - 1);
        for (size_t i = 0; i < test_size; i++) {
            size_t label = label_dist(generator);
            labels.emplace_back(label);
        }
        return labels;
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


    template<int dim>
    auto randomPoints(size_t test_size) -> std::vector<pyscan::Point<dim>> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0, 1.0);

        std::vector<pyscan::Point<dim>> points(test_size, pyscan::Point<dim>());
        for (size_t i = 0; i < test_size; i ++) {
            for (size_t j = 0; j < dim + 1; j++) {
                points[i][j] = distribution(generator);
            }
            points[i][dim] = 1.0;
        }
        return points;
    }

    template<int dim>
    auto randomWPoints(size_t test_size) -> std::vector<pyscan::WPoint<dim>> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);

        std::vector<pyscan::WPoint<dim>> points(test_size, pyscan::WPoint<dim>());
        for (size_t i = 0; i < test_size; i ++) {
            for (size_t j = 0; j < dim; j++) {
                points[i][j] = distribution(generator);
            }
            points[i][dim] = 1.0;
            points[i].set_weight(1.0);
        }
        return points;
    }

    template<int dim>
    auto randomLPoints(size_t test_size, size_t label_count) -> std::vector<pyscan::LPoint<dim>> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);

        auto labels = randomLabels(test_size, label_count);
        std::vector<pyscan::LPoint<dim>> points(test_size, pyscan::LPoint<dim>());

        for (size_t i = 0; i < test_size; i ++) {
            for (size_t j = 0; j < dim; j++) {
                points[i][j] = distribution(generator);
            }
            points[i][dim] = 1.0;
            points[i].set_weight(1.0);
            points[i].set_label(labels[i]);
        }
        return points;
    }

    auto randomPoints2(size_t test_size) -> std::vector<pyscan::Point<2>> {
        return randomPoints<2>(test_size);
    }

    auto randomPoints3(size_t test_size) -> std::vector<pyscan::Point<3>> {
        return randomPoints<3>(test_size);
    }
    auto randomWPoints2(size_t test_size) -> std::vector<pyscan::WPoint<2>> {
        return randomWPoints<2>(test_size);
    }

    auto randomWPoints3(size_t test_size) -> std::vector<pyscan::WPoint<3>> {
        return randomWPoints<3>(test_size);
    }


    auto randomLPoints2(size_t test_size, size_t label_count) -> std::vector<pyscan::LPoint<2>> {
        return randomLPoints<2>(test_size, label_count);
    }

    auto randomLPoints3(size_t test_size, size_t label_count) -> std::vector<pyscan::LPoint<3>> {
        return randomLPoints<3>(test_size, label_count);
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


    template < int dim>
    auto randomLPointsUnique(size_t test_size) -> std::vector<pyscan::LPoint<dim>> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);

        std::vector<pyscan::LPoint<dim>> points(test_size, pyscan::LPoint<dim>());

        for (size_t i = 0; i < test_size; i ++) {
            for (size_t j = 0; j < dim; j++) {
                points[i][j] = distribution(generator);
            }
            points[i][dim] = 1.0;
            points[i].set_weight(1.0);
            points[i].set_label(i);
        }
        return points;
    }

    auto randomLPointsUnique2(size_t test_size) -> std::vector<pyscan::LPoint<>> {
        return randomLPointsUnique<2>(test_size);
    }

    auto randomLPointsUnique3(size_t test_size) -> std::vector<pyscan::LPoint<3>> {
        return randomLPointsUnique<3>(test_size);
    }


    auto removeLabels(pyscan::lpoint_list_t const& pts) -> pyscan::wpoint_list_t {

        pyscan::wpoint_list_t wlist(pts.size(), pyscan::WPoint<>());
        std::transform(pts.begin(), pts.end(), wlist.begin(), [&](pyscan::LPoint<> const& lpt){
            return pyscan::WPoint<>(lpt.get_weight(), lpt[0], lpt[1], lpt[2]);
        });
        return wlist;
    }

    auto removeLW(pyscan::lpoint_list_t const& pts) -> pyscan::point_list_t {

        pyscan::point_list_t plist(pts.size(), pyscan::Point<>());
        std::transform(pts.begin(), pts.end(), plist.begin(), [&](pyscan::LPoint<> const& lpt){
            return pyscan::Point<>(lpt[0], lpt[1], lpt[2]);
        });
        return plist;
    }

    auto addWeights(pyscan::point_list_t const& pts) -> pyscan::wpoint_list_t {
        pyscan::wpoint_list_t list;
        for_each(pts.begin(), pts.end(), [&](auto const& pt){
            list.emplace_back(pyscan::WPoint<>(1.0, pt[0], pt[1], pt[2]));
        });
        return list;
    }

    auto addLabels(pyscan::wpoint_list_t const& pts, pyscan::label_list_t const& labels) -> pyscan::lpoint_list_t {
        pyscan::lpoint_list_t list;
        auto l_it = labels.begin();
        for_each(pts.begin(), pts.end(), [&](auto const& pt){
            list.emplace_back(pyscan::LPoint<>(1.0, pt.get_weight(), pt[0], pt[1], pt[2]));
            l_it++;
        });
        return list;
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
