//
// Created by mmath on 12/10/18.
//

#ifndef PYSCAN_SAMPLING_HPP
#define PYSCAN_SAMPLING_HPP

#include <random>
#include <vector>
#include <unordered_map>
#include <functional>
#include <queue>


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

    template<class It, class URNG, class OutIt>
    OutIt random_sample_wor(It b, It e,
                           URNG &&g,
                           const std::function<double(decltype(*b))> &wf,
                           size_t sample_size, OutIt out) {
        using el_t = decltype(*b);
        using node_t = std::tuple<double, el_t>;

        auto cmp = [](const node_t &e1, const node_t &e2) {
            return std::get<0>(e1) < std::get<0>(e2);
        };

        std::uniform_real_distribution<double> dist(0, 1);
        std::priority_queue<node_t, std::vector<node_t>, decltype(cmp)> sample_heap(cmp);
        for (; b != e; ++b) {

            sample_heap.push(std::make_tuple(-std::pow(dist(g), 1.0 / wf(*b)), *b));
            if (sample_heap.size() >= sample_size) {
                //Remove smallest element
                sample_heap.pop();
            }
        }

        while (!sample_heap.empty()) {
            *out = std::get<1>(sample_heap.top());
            ++out;
            sample_heap.pop();
        }
        return out;
    }


    template<class It, class URNG, class OutIt>
    OutIt random_sample_wr(It b, It e,
                            URNG &&g,
                            std::function<double(It)> wf,
                            size_t sample_size, OutIt out) {

        std::vector<double> cum_sums;
        double initial = 0;
        for (It bb = b; bb != e; bb++) {
            initial += wf(bb);
            cum_sums.emplace_back(initial);
        }

        std::vector<size_t> indices;
        std::uniform_real_distribution<double> dist(0, initial);
        for (size_t s = 0; s < sample_size; s++) {

            auto it = std::upper_bound(cum_sums.begin(), cum_sums.end(), dist(g));
            if (it == cum_sums.end()) {
                it--;
            }
            indices.emplace_back(it - cum_sums.begin());
        }
        std::sort(indices.begin(), indices.end());
        auto ib = indices.begin();
        size_t i = 0;
        while(ib != indices.end()) {
            while (i != *ib) {
                i++;
                b++;
            }
            *out = b;
            ib++;
            out++;
        }
        return out;
    }

    template<class It, class URNG>
    auto random_w_selection(It b, It e, URNG&& gen, const std::function<double(decltype(*b))> &wf) -> decltype(*b) {
        std::vector<double> cum_sums;
        cum_sums.reserve(e - b);
        double initial = 0;
        for (auto bb = b; bb != e; bb++) {
            initial += wf(*bb);
            cum_sums.emplace_back(initial);
        }

        std::uniform_real_distribution<double> dist(0, initial);
        auto it = std::upper_bound(cum_sums.begin(), cum_sums.end(), dist(gen));
        if (it == cum_sums.end()) {
            it--;
        }
        return  *(b + (it - cum_sums.begin()));
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

    template<typename It, typename Cmp>
    It weighted_median_slow(It b, It e, Cmp cmp, const std::function<double(decltype(*b))>& wf) {
        auto total_weight = [&wf](It b_local, It e_local) {
            double tw = 0;
            for (auto b1 = b_local; b1 != e_local; ++b1) {
                tw += wf(*b1);
            }
            return tw;
        };

        auto b_store = b;
        std::sort(b, e, cmp);
        double tw = total_weight(b, e);
        double tmp_w = 0;
        for (; b != e; ++b) {
            tmp_w += wf(*b);
            if (tmp_w >= tw / 2.0) {
                break;
            }
        }

        if (b_store == b) {
            return b;
        } else if ( tw / 2.0 -  (tmp_w - wf(*(b))) < tmp_w - tw / 2.0) {
            return b - 1;
        } else {
            return b;
        }
    }

    template<typename It, typename Cmp, typename URNG>
    It weighted_median(It b, It e, Cmp cmp, const std::function<double(decltype(*b))>& wf, URNG &&gen) {
        auto total_weight = [&wf](It b_local, It e_local) {
            double tw = 0;
            for (auto b1 = b_local; b1 != e_local; ++b1) {
                tw += wf(*b1);
            }
            return tw;
        };
        double lw = 0;
        double rw = 0;

        size_t min_buff_size = 10;
        while (e - b > min_buff_size) {
            std::vector<decltype(*b)> output_buffer;
            output_buffer.reserve(10);

            auto it = random_sample_wor(b, e, gen, 10, wf, std::back_inserter(output_buffer));
            auto split_pt = weighted_median_slow(output_buffer.begin(), it, cmp, wf);
            auto splt_pt = std::partition(b, e, [&](decltype(*it) const& p) {
                return cmp(p, *split_pt);
            });

            double tmp_lw = total_weight(b, split_pt) + lw;
            double tmp_rw = total_weight(split_pt, e) + rw;
            if (tmp_lw < tmp_rw) {
                e = split_pt;
                lw = tmp_lw;
            } else {
                b = split_pt;
                rw = tmp_rw;
            }
        }
        return weighted_median_slow(b, e, cmp, wf);
    }



}
#endif //PYSCAN_SAMPLING_HPP
