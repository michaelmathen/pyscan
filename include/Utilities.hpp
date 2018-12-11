//
// Created by mmath on 6/28/17.
//

#ifndef PYSCAN_UTILITIES_HPP
#define PYSCAN_UTILITIES_HPP
#include <vector>
#include <algorithm>
#include <sstream>
#include <ostream>
#include <cstdint>
#include <limits>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <random>
#include <numeric>

namespace pyscan {


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

    namespace util {

        const int64_t max_float_distance = 5;
        const int64_t max_int64 = std::numeric_limits<int64_t>::max();

        inline int64_t ulpsDistance(const double a, const double b) {
            if (a == b) return 0;

            if (std::isnan(a) || std::isnan(b)) return max_int64;
            if (std::isinf(a) || std::isinf(b)) return max_int64;

            static_assert(sizeof(double) == sizeof(int64_t));
            int64_t ia, ib;
            memcpy(&ia, &a, sizeof(double));
            memcpy(&ib, &b, sizeof(double));

            if ((ia < 0) != (ib < 0)) return max_int64;

            return std::abs(ia - ib);
        }

        
        inline bool aeq(double a, double b) {
            //return std::abs(a - b) / std::max(std::abs(a), std::abs(b)) <= .0000000001;
            return ulpsDistance(a, b) < max_float_distance;
        }

        inline bool alt(double a, double b) {
            return a < b && !aeq(a, b);
        }

        inline bool alte(double a, double b) {
            return a <= b || aeq(a, b);
        }

        inline double det2(double a1, double b1, double a2, double b2) {
            return a1 * b2 - a2 * b1;
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
    }

}
#endif //PYSCAN_UTILITIES_HPP
