//
// Created by Michael on 7/12/2018.
//

#ifndef PYSCAN_SPARSEGRID_H
#define PYSCAN_SPARSEGRID_H

#include <cstdint>
#include <vector>

#include "Point.hpp"
#include "Utilities.hpp"

namespace pyscan {

    inline uint64_t to_code(Point<2> const& pt, uint32_t r, double scale, double min_x, double min_y) {
        double x = getX(pt), y = getY(pt);
        auto a = static_cast<uint32_t>((x - min_x) / scale * r),
                b = static_cast<uint32_t>((y - min_y) / scale * r);
        return morton(a, b);
    }

    template<typename T>
    class SparseGrid {
        std::vector<uint64_t> z_vals;
        std::vector<T> z_order_pts;

        double mnx;
        double mny;
        double scale;
        uint32_t r;



    public:
        using pt_it = decltype(z_order_pts.begin());
        using coord_func_t = std::function<double(T)>;


        SparseGrid(std::vector<T> const& items, int32_t r) : z_order_pts(items), r(r) {


            pt_it min_x, max_x;
            std::tie(min_x, max_x) = std::minmax_element(z_order_pts.begin(),
                                                         z_order_pts.end(), [&](T const& p1, T const& p2) {
                        return getX(p1) < getX(p2);
                    });


            pt_it min_y, max_y;

            std::tie(min_y, max_y) = std::minmax_element(z_order_pts.begin(), z_order_pts.end(),
                                                         [&](T const& p1, T const& p2) {
                                                             return getY(p1) < getY(p2);
                                                         });

            scale = std::max(getX(*max_x) - getX(*min_x), getY(*max_y) - getY(*min_y));
            mnx = getX(*min_x);
            mny = getY(*min_y);

            z_vals.reserve(z_order_pts.size());
            for (auto& pt : z_order_pts) {
                z_vals.push_back(to_code(pt, r, scale, mnx, mny));
            }
            auto perm = sort_permutation(z_vals.begin(), z_vals.end(), std::less<uint32_t>());
            apply_permutation_in_place(z_order_pts, perm);
            apply_permutation_in_place(z_vals, perm);
        }

        auto operator()(int32_t i, int32_t j) -> std::tuple<pt_it, pt_it> {
            uint64_t code = morton(i, j);
            auto begin = std::lower_bound(z_vals.begin(), z_vals.end(), code);
            auto end = std::upper_bound(begin, z_vals.end(), code);
            return std::make_tuple(z_order_pts.begin() + (begin - z_vals.begin()), z_order_pts.begin() + (end - z_vals.begin()));
        }

        auto get_resolution() -> double {
            return 1 / static_cast<double>(r);
        }
    };
}

#endif //PYSCAN_SPARSEGRID_H
