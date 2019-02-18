//
// Created by mmath on 12/10/18.
//

#ifndef PYSCAN_SAMPLING_HPP
#define PYSCAN_SAMPLING_HPP

#include <random>
#include <vector>
#include <unordered_map>
#include <functional>


namespace pyscan {
    template<class Vect, class URNG>
    auto random_sample_wor(Vect const &arr,
                           URNG &&g,
                           const std::function<double(decltype(*arr.begin()))> &wf,
                           size_t sample_size) -> Vect {
        using el_t = typename Vect::value_type;
        using node_t = std::tuple<double, el_t>;

        auto cmp = [](const node_t &e1, const node_t &e2) {
            return std::get<0>(e1) < std::get<0>(e2);
        };

        std::uniform_real_distribution<double> dist(0, 1);
        std::priority_queue<node_t, std::vector<node_t>, decltype(cmp)> sample_heap(cmp);
        for (size_t i = 0; i < arr.size(); ++i) {

            sample_heap.push(std::make_tuple(-std::pow(dist(g), 1.0 / wf(arr[i])), arr[i]));
            if (sample_heap.size() >= sample_size) {
                //Remove smallest element
                sample_heap.pop();
            }
        }

        Vect output;
        while (!sample_heap.empty()) {
            output.emplace_back(std::get<1>(sample_heap.top()));
            sample_heap.pop();
        }
        return output;
    }

    template<class Vect, class URNG>
    auto random_sample_wr(Vect const &arr, URNG &&g, size_t sample_size) -> Vect {
        Vect output;
        std::uniform_int_distribution<size_t> d(0, arr.size() - 1);
        for (size_t i = 0; i < sample_size; ++i) {
            output.push_back(arr[d(g)]);
        }
        return output;
    }

    template<class Vect, class URNG>
    auto random_sample_wor(Vect const &arr, URNG &&g, size_t sample_size) -> Vect {
        Vect output;
        std::unordered_map<size_t, decltype(*arr.begin())> key;
        sample_size = std::min(arr.size(), sample_size);
        for (size_t i = 0; i < sample_size; ++i) {
            std::uniform_int_distribution<decltype(i)> d(i, arr.size() - 1);
            auto rix = d(g);
            decltype(*arr.begin()) rel;
            if (key.find(rix) == key.end()) {
                rel = arr[rix];
            } else {
                rel = key[rix];
            }
            decltype(*arr.begin()) arr_el;
            if (key.find(i) == key.end()) {
                arr_el = arr[i];
            } else {
                arr_el = key[i];
                key.erase(i); // no chance we will access this again.
            }
            output.push_back(rel);
            key[rix] = arr_el;
        }
        return output;
    }

    template<class Vect, class URNG>
    Vect random_sample_wor(Vect &arr, URNG &&g, size_t sample_size) {
        sample_size = std::min(arr.size(), sample_size);
        for (size_t i = 0; i < sample_size; ++i) {
            std::uniform_int_distribution<decltype(i)> d(i, arr.size() - 1);
            std::swap(arr[i], arr[d(g)]);
        }
        return Vect(arr.begin(), arr.begin() + sample_size);
    }
}
#endif //PYSCAN_SAMPLING_HPP
