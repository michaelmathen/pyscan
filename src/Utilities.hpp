//
// Created by mmath on 6/28/17.
//

#ifndef PYSCAN_UTILITIES_HPP
#define PYSCAN_UTILITIES_HPP
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <ostream>
#include <cstdint>

namespace pyscan {


    double invsqrt(double number);


    template<class T>
    auto operator<<(std::ostream &os, const T &t) -> decltype(t.print(os), os) {
        t.print(os);
        return os;
    }


    template<typename T>
    std::ostream &operator<<(std::ostream &out, std::vector<T> const &els) {
        out << "[";
        for (auto i : els) {
            out << i << ", ";
        }
        out << "]";
        return out;
    }

    template<typename T, typename Compare>
    std::vector<std::size_t> sort_permutation(const std::vector<T> &vec, Compare compare) {
        std::vector<std::size_t> p(vec.size());
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(),
                  [&](std::size_t i, std::size_t j) { return compare(vec[i], vec[j]); });
        return p;
    }

    template<typename It, typename Compare>
    std::vector<std::size_t> sort_permutation(It begin, It end, Compare compare) {
        std::vector<std::size_t> p(end - begin);
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) {
            return compare(*(begin + i), *(begin + j));
        });
        return p;
    }

    template<typename T>
    void apply_permutation_in_place(std::vector<T> &vec, std::vector<std::size_t> const &p) {
        std::vector<bool> done(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i) {
            if (done[i]) {
                continue;
            }
            done[i] = true;
            std::size_t prev_j = i;
            std::size_t j = p[i];
            while (i != j) {
                std::swap(vec[prev_j], vec[j]);
                done[j] = true;
                prev_j = j;
                j = p[j];
            }
        }
    }


    template<typename P_it, typename I_it, typename W_it, typename Comp>
    void quantiles(P_it begin, P_it end, I_it i_begin, I_it i_end, W_it w_begin, Comp comp) {
        auto order = sort_permutation(begin, end, comp);
        int r = i_end - i_begin;
        double total_weight = std::accumulate(w_begin, w_begin + (end - begin), 0, std::plus<>());

        double eps_s = total_weight / r;
        double curr_weight = 0;
        std::for_each(order.begin(), order.end(), [&](size_t ix) {
            curr_weight += *(w_begin + ix);
            if (curr_weight > eps_s) {
                *i_begin = *(ix + begin);
                i_begin++;
                curr_weight = 0;
            }
        });
    }



    uint64_t morton(uint32_t a, uint32_t b) {

        uint64_t x = a, y = b;
        x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
        x = (x | (x << 8)) & 0x00FF00FF00FF00FF;
        x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0F;
        x = (x | (x << 2)) & 0x3333333333333333;
        x = (x | (x << 1)) & 0x5555555555555555;

        y = (y | (y << 16)) & 0x0000FFFF0000FFFF;
        y = (y | (y << 8)) & 0x00FF00FF00FF00FF;
        y = (y | (y << 4)) & 0x0F0F0F0F0F0F0F0F;
        y = (y | (y << 2)) & 0x3333333333333333;
        y = (y | (y << 1)) & 0x5555555555555555;

        return x | (y << 1);
    }

}
#endif //PYSCAN_UTILITIES_HPP
