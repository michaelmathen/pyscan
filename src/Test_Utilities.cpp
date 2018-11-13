//
// Created by mmath on 10/2/17.
//

#include <random>
#include <algorithm>
#include <Point.hpp>


#include "Test_Utilities.hpp"

namespace pyscantest {

    auto randomPoints(int test_size) -> std::vector<pyscan::pt2_t> {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution (0.0,1.0);
        std::vector<pyscan::pt2_t> points;
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


    auto randomLabels(int test_size, size_t num_labels) -> pyscan::label_list_t {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::vector<size_t> labels;
        using label_t = size_t;

        std::uniform_int_distribution<label_t> label_dist(0, num_labels - 1);
        for (int i = 0; i < test_size; i++) {
            size_t label = label_dist(generator);
            labels.emplace_back(label);
        }
        return labels;
    }


    auto randomLPoints(int test_size, int label_count) -> std::vector<pyscan::lpt2_t> {
        auto pts = randomVec(test_size);
        auto lbls = randomLabels(test_size, label_count);
        std::vector<pyscan::lpt2_t> lpoints(test_size, pyscan::lpt2_t ());
        auto wp = lpoints.begin();
        auto lp = lbls.begin();
        std::for_each(pts.begin(), pts.end(), [&](auto const& pt) {
            *wp = pyscan::LPoint<>(*lp, 1.0, pt[0], pt[1], 1.0);
            lp++;
            wp++;
        });
        return lpoints;
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



    auto randomLPointsUnique(int test_size) -> std::vector<pyscan::LPoint<>> {
        auto pts = randomVec(test_size);
        std::vector<pyscan::LPoint<>> lpoints(test_size, pyscan::LPoint<>());
        auto wp = lpoints.begin();
        size_t curr_label = 0;
        std::for_each(pts.begin(), pts.end(), [&](auto const& pt) {
            *wp = pyscan::LPoint<>(curr_label, 1.0, pt[0], pt[1], 1.0);
            curr_label++;
            wp++;
        });
        return lpoints;
    }

    auto removeLabels(pyscan::lpoint_list_t const& pts) -> pyscan::wpoint_list_t {

        pyscan::wpoint_list_t wlist(pts.size(), pyscan::WPoint<>());
        std::transform(pts.begin(), pts.end(), wlist.begin(), [&](pyscan::LPoint<> const& lpt){
            return pyscan::WPoint<>(lpt.get_weight(), lpt[0], lpt[1], lpt[2]);
        });
        return wlist;
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
