#ifndef __SPARSE_GRID_H__
#define __SPARSE_GRID_H__

#include <iostream>
#include <assert.h>
#include <map>
#include <algorithm>
#include <vector>
#include <type_traits>
#include <cmath>
#include <unordered_map>

/*
 * This class implements high level behaviour of a regular grid over R^2. We can use this to implement features of
 * other grids needed for this code.
 */
template <typename Derived, typename Storage_t>
class Grid {

public:
    using value_type = typename Storage_t::value_type;
    using Iterator_t = typename Storage_t::const_iterator;

    Grid(double bx, double by, double sc, double m_r) :
        bx(bx),
        by(by),
        scale(sc),
        min_res(m_r),
        r(static_cast<uint32_t>(lround(floor(sc / m_r)))) {

    }

    double get_resolution() const {
        return scale / static_cast<double>(r);
    }

    std::pair<uint32_t, uint32_t> get_cell(double x, double y) const {
        if (x - bx < 0) {
            std::cout << x << " is less than " << bx << std::endl;
        }
        if (y - by < 0){
            std::cout << y << " is less than " << by << std::endl;
        }
        assert(x - bx >= 0);
        assert(y - by >= 0);
        auto a = static_cast<uint32_t>((x - bx) / scale * r);
        auto b = static_cast<uint32_t>((y - by) / scale * r);
        a = (a == r) ? r - 1 : a;
        b = (b == r) ? r - 1 : b;
        return std::make_pair(a, b);
    }

    std::tuple<double, double> get_center(uint32_t ix, uint32_t iy) const {
        return std::make_tuple((2.0 * ix + 1) * min_res /2  + bx, (2.0 * iy + 1) * min_res /2  + by);
    }

    std::pair<uint32_t, uint32_t> get_cell_offset(double x, double y) const {
        /*
         * This method ignores the lower bound offset of the grid, so that we
         * can just use the offset into the grid.
         */
        auto a = static_cast<uint32_t>(x / scale * r);
        auto b = static_cast<uint32_t>(y / scale * r);
        a = (a == r) ? r - 1 : a;
        b = (b == r) ? r - 1 : b;
        return std::make_pair(a, b);
    }

    std::tuple<double, double> get_lower_corner(uint32_t ix, uint32_t iy) const {
        return std::make_tuple(ix * min_res + bx, iy * min_res + by);
    }

    std::tuple<double, double> get_lower_corner(uint64_t code) const {
        if (r == 0)
            return std::make_tuple(0, 0);
        uint32_t x_i = code % r;
        uint32_t y_i = code / r;
        return get_lower_corner(x_i, y_i);
    }

    std::tuple<double, double> get_lower_corner() const {
        return std::make_tuple(bx, by);
    }

    std::tuple<double, double> get_center(uint64_t code) const {
        if (r == 0)
            return std::make_tuple(0, 0);
        uint32_t x_i = code % r;
        uint32_t y_i = code / r;
        return get_center(x_i, y_i);
    }

    template <typename V>
    uint64_t get_code(const V& pt) const {
        auto res = get_cell(pt);
        return res.second * r + res.first;
    }

    uint64_t get_code(double x, double y) const {
        auto res = get_cell(x, y);
        return res.second * r + res.first;
    }

    template <typename V>
    std::pair<uint32_t, uint32_t> get_cell(const V& pt) const {
        return get_cell(pt(0), pt(1));
    }

    std::tuple<uint32_t, uint32_t> get_cell(uint64_t code) const {
        uint32_t x_i = code % r;
        uint32_t y_i = code / r;
        return std::make_tuple(x_i, y_i);
    }

    Iterator_t begin() const {
        return derived_const().implement_begin();
    }

    Iterator_t end() const {
        return derived_const().implement_end();
    }

    uint32_t get_grid_size() const {
        return r;
    }

protected:
    //The bottom x and bottom y corner of this grid.
    double bx, by;
    //The scale that maps points in doubles to points in the grid.
    double scale;
    double min_res;
    uint32_t r;
private:
    Derived const& derived_const() const {
        return *static_cast<Derived const*>(this);
    }

    Derived& derived() {
        return *static_cast<Derived*>(this);
    }

};

template <typename T>
class SparseGrid : public Grid<SparseGrid<T>, std::multimap<uint64_t, T>> {
public:
    using bbox_t = std::tuple<double, double, double, double>;

    using iterator_t = typename std::multimap<uint64_t, T>::const_iterator;


    SparseGrid(bbox_t bb, const std::vector<T>& items, double min_res) :
        Grid<SparseGrid<T>, std::multimap<uint64_t, T>>(std::get<0>(bb), std::get<1>(bb),
                std::max(std::abs(std::get<0>(bb) - std::get<2>(bb)),
                       std::abs(std::get<1>(bb) - std::get<3>(bb))),
                       min_res) {
        for (auto& pt: items) {
            z_pts.emplace(this->get_code(pt), pt);
        }
    }

    SparseGrid(bbox_t bb, double min_res) :
        Grid<SparseGrid<T>, std::multimap<uint64_t, T>>(std::get<0>(bb), std::get<1>(bb),
                std::max(std::abs(std::get<0>(bb) - std::get<2>(bb)),
                       std::abs(std::get<1>(bb) - std::get<3>(bb))),
                       min_res) {
    }

    std::pair<iterator_t, iterator_t> operator()(uint32_t i, uint32_t j) const {
        if ((i >= this->r) || (j >= this->r)) {
            return std::make_pair(z_pts.end(), z_pts.end());
        }
        uint64_t code = j * this->r + i;
        return z_pts.equal_range(code);
    }

    std::pair<iterator_t, iterator_t> operator()(uint64_t code) const {
        return z_pts.equal_range(code);
    }

    iterator_t implement_begin() const {
        return z_pts.cbegin();
    }
    iterator_t implement_end() const {
        return z_pts.cend();
    }
    std::multimap<uint64_t, T> z_pts;
};


class WeightedGrid : public Grid<WeightedGrid, std::unordered_map<uint64_t, double>> {
public:
    using bbox_t = std::tuple<double, double, double, double>;

    template<typename It, typename Wf, typename Xf>
    WeightedGrid(bbox_t bb, It b, It e, double min_res, Wf wf, Xf xf) :
            Grid<WeightedGrid, std::unordered_map<uint64_t, double> >(std::get<0>(bb), std::get<1>(bb),
                                std::max(std::abs(std::get<0>(bb) - std::get<2>(bb)),
                                         std::abs(std::get<1>(bb) - std::get<3>(bb))),
                                min_res) {
        assert(this->r > 0);
        for (auto it = b; it != e; it++) {
            auto [ix, iy] = xf(*it);
            if (weights.find(get_code(ix, iy)) == weights.end()) {
                weights.emplace(get_code(ix, iy), wf(*it));
            }else {
                weights[get_code(ix, iy)] += wf(*it);
            }
        }
    }


    double operator()(uint32_t i, uint32_t j) const {
        uint64_t code = j * r + i;
        auto it = weights.find(code);
        if (it == weights.end()) return 0.0;
        return it->second;
    }

    Iterator_t implement_begin() const {
        return weights.cbegin();
    }

    Iterator_t implement_end() const {
        return weights.cend();
    }

    std::unordered_map<uint64_t, double> weights;
};


using bbox_t = std::tuple<double, double, double, double>;

/*
 * Invokes the scan function on every corner point of a grid.
 */
template<typename Scan_f>
void grid_scan(bbox_t bbox, double grid_res, Scan_f func) {
    auto [mnx, mny, mxx, mxy] = bbox;
    for (double x = mnx; x < mxx; ) {
        for (double y = mny; y < mxy;) {
            func(x, y);
            y = y + grid_res;
        }
        x = x + grid_res;
    }
}


/*
 * Invokes the scan function on a subset of points.
 */
#endif
