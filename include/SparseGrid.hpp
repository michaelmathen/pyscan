#ifndef __SPARSE_GRID_H__
#define __SPARSE_GRID_H__

#include <assert.h>
#include <map>
#include <algorithm>
#include <vector>

template <typename T>
class SparseGrid {
public:
    using iterator_t = typename std::multimap<uint64_t, T>::const_iterator;

    SparseGrid(const std::vector<T>& items, uint32_t grid_r) : r(grid_r) {
        assert(r > 0);
        if (items.size() == 0) {
            return;
        }
        
        for (auto& pt: items) {
            z_pts.emplace(get_code(pt), pt);
        }
    }

    std::pair<iterator_t, iterator_t> operator()(uint32_t i, uint32_t j) const {
        if ((i >= r) || (j >= r)) {
            return std::make_pair(z_pts.end(), z_pts.end());
        }
        uint64_t code = j * r + i;
        return z_pts.equal_range(code);
    }

    double get_resolution() const {
        return 1.0 / static_cast<double>(r);
    }

    std::pair<uint32_t, uint32_t> get_cell(const T& pt) const {
        uint32_t a = static_cast<uint32_t>(pt(0) * r);
        uint32_t b = static_cast<uint32_t>(pt(1) * r);
        a = (a == r) ? r - 1 : a;
        b = (b == r) ? r - 1 : b;
        return std::make_pair(a, b);
    }

    uint64_t get_code(const T& pt) const {
        auto res = get_cell(pt);
        return res.second * r + res.first;
    }

    iterator_t begin() const {
        return z_pts.cbegin();
    }

    iterator_t end() const {
        return z_pts.cend();
    }
    
private:
    std::multimap<uint64_t, T> z_pts;
    uint32_t r;
};

#endif
