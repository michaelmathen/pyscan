//
// Created by Michael on 7/12/2018.
//

#ifndef PYSCAN_SPARSEGRID_H
#define PYSCAN_SPARSEGRID_H

#include <cstdint>
#include <vector>
#include <iostream>
#include <map>

#include "Point.hpp"
#include "Utilities.hpp"

namespace pyscan {

    inline std::tuple<int32_t, int32_t> to_cell(Point<2> const& pt, int32_t r, double scale, double min_x, double min_y) {
        double x = getX(pt),
                y = getY(pt);
        auto a = static_cast<int32_t>((x - min_x) / scale * r),
                b = static_cast<int32_t>((y - min_y) / scale * r);
        a = a == r ? r - 1 : a;
        b = b == r ? r - 1 : b;
        return std::make_tuple(a, b);
    }

    inline int64_t to_code(Point<2> const& pt, int32_t r, double scale, double min_x, double min_y) {
        int64_t a, b;
        std::tie(a, b) = to_cell(pt, r, scale, min_x, min_y);
        return b * r + a;
    }

    template<typename T>
    class SparseGrid {
        std::multimap<int64_t, T> z_pts;
        double mnx;
        double mny;
        double scale;
        int32_t r;

    public:

        using buck_it = decltype(z_pts.begin());

        SparseGrid(std::vector<T> const& items, int32_t grid_r) : r(grid_r) {
            if (items.size() == 0) {
                return;
            }

            using pt_it = decltype(items.begin());
            pt_it min_x, max_x;
            std::tie(min_x, max_x) = std::minmax_element(items.begin(), items.end(), [&](T const& p1, T const& p2) {
                        return getX(p1) < getX(p2);
            });
            pt_it min_y, max_y;
            std::tie(min_y, max_y) = std::minmax_element(items.begin(), items.end(),
                                                         [&](T const& p1, T const& p2) {
                                                             return getY(p1) < getY(p2);
                                                         });
            scale = std::max(getX(*max_x) - getX(*min_x), getY(*max_y) - getY(*min_y));
            mnx = getX(*min_x);
            mny = getY(*min_y);
            for (auto& pt : items) {
                z_pts.emplace(to_code(pt, r, scale, mnx, mny), pt);
            }
        }

        auto operator()(int32_t i, int32_t j) -> std::tuple<buck_it, buck_it> {

            if ((i < 0) || (i > r) || (j < 0) || (j > r)) {
                return std::make_tuple(z_pts.end(), z_pts.end());
            }
            int64_t code = i * r + j;
            return z_pts.equal_range(code);
        }

        auto get_resolution() -> double {
            return 1 / static_cast<double>(r);
        }


        std::tuple<int32_t, int32_t> get_cell(T const& pt) {
            return to_cell(pt, r, scale, mnx, mny);
        }

        auto begin() -> decltype(z_pts.begin()) {
            return z_pts.begin();
        }

        auto end() -> decltype(z_pts.begin()) {
            return z_pts.end();
        }

    };
}

#endif //PYSCAN_SPARSEGRID_H
