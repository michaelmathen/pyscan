//
// Created by mmath on 5/28/17.
//

#ifndef PYSCAN_RECTANGLESCAN_HPP
#define PYSCAN_RECTANGLESCAN_HPP
#include <memory>
#include <vector>
#include <forward_list>
#include <algorithm>
#include <cmath>
#include <random>
#include <deque>
#include <functional>
#include <tuple>
#include <limits>
#include <iostream>
#include "BloomFilter.hpp"

#include "Range.hpp"
#include "Statistics.hpp"
#include "Point.hpp"
#include "Utilities.hpp"

//#define DEBUG

namespace pyscan {
    using point_it= std::vector<Point<> >::iterator;
    using subgrid = std::tuple<int, int, int, int, double>;


    class Subgrid {
        size_t u_x, u_y, l_x, l_y;
        double value;
    public:
        Subgrid(size_t ux, size_t uy, size_t lx, size_t ly, double val) :
                u_x(ux),
                u_y(uy),
                l_x(lx),
                l_y(ly),
                value(val) {}


        std::string toString() const {
           std::stringstream ss;
           ss << "(" << u_x << " " << u_y << " " << l_x << " " << l_y << " " << value << ")";
           return ss.str();
        }

        void print(std::ostream& os) const {
           os <<  "(" << u_x << " " << u_y << " " << l_x << " " << l_y << " " << value << ")";
        }

        bool contains(Point<> const& pt) const {
            return u_x >= pt(0) && pt(0) >= l_x && u_y >= pt(1) && pt(1) >= l_y;
        }

        void setValue(double v) {
          value = v;
        }
        size_t lowX() const { return l_x; }
        size_t upX() const { return u_x; }
        size_t lowY() const { return l_y; }
        size_t upY() const { return u_y; }
        double fValue() const { return value; }
    };



    class Rectangle : public Range<2> {
        double u_x, u_y, l_x, l_y;
    public:
        Rectangle(const pt2_t& p1, const pt2_t& p2, const pt2_t& p3, const pt2_t& p4) {
            u_x = std::max({p1(0), p2(0), p3(0), p4(0)});
            l_x = std::min({p1(0), p2(0), p3(0), p4(0)});
            u_y = std::max({p1(1), p2(1), p3(1), p4(1)});
            l_y = std::min({p1(1), p2(1), p3(1), p4(1)});
        }

        Rectangle() : u_x(0.0), u_y(0.0), l_x(0.0), l_y(0.0) {}

        Rectangle(double ux, double uy, double lx, double ly) : u_x(ux), u_y(uy), l_x(lx), l_y(ly) {}

        inline bool contains(const pt2_t& p1) const final {
            return (u_x > p1(0) && p1(0) >= l_x && u_y > p1(1) && p1(1) >= l_y);
        }

        inline bool intersects_segment(const pt2_t &p1, const pt2_t &p2) const final {
            //Check if any of points are inside of the rectangle.
            if (contains(p1) || contains(p2)) {
                return true;
            } else {
                //Check if any of the segments cross the line segment defined by p1<->p2.
                auto ur_pt = pt2_t(u_x, u_y, 1.0);
                auto ul_pt = pt2_t(l_x, u_y, 1.0);
                auto lr_pt = pt2_t(u_x, l_y, 1.0);
                auto ll_pt = pt2_t(l_x, l_y, 1.0);
                return crosses_segment(ur_pt, ul_pt, p1, p2) || crosses_segment(ul_pt, ll_pt, p1, p2) ||
                        crosses_segment(ur_pt, lr_pt, p1, p2) || crosses_segment(lr_pt, ll_pt, p1, p2);
            }
        }

        std::string toString() const {
            std::stringstream ss;
            ss << "Rectangle(" << u_x << ", " << u_y << ", " << l_x << ", " << l_y << ")";
            return ss.str();
        }
        double lowX() const { return l_x; }
        double upX() const { return u_x; }
        double lowY() const { return l_y; }
        double upY() const { return u_y; }
    };

    enum class I_Type {
        VALUE,
        MAX
    };

    class Grid {
        /*
         * Grid as defined in the SODA paper. Construction takes O(mlog r + r^2) time where
         * m = end - begin.
         */
        long r;
        std::vector<double> red_counts;
        std::vector<double> blue_counts;
        std::vector<double> x_coords;
        std::vector<double> y_coords;
        double total_red_weight = 0;
        double total_blue_weight = 0;
    public:
        Grid(size_t r_arg, wpoint_list_t const& red, wpoint_list_t const& blue);
        Grid(wpoint_list_t const& red_points, wpoint_list_t const& blue_points);
        Grid(point_list_t const& net, wpoint_list_t const& red, wpoint_list_t const& blue);
        double totalRedWeight() const;
        double totalBlueWeight() const;
        double redCount(size_t row, size_t col) const;
        double blueCount(size_t row, size_t col) const;
        double redWeight(size_t row, size_t col) const;
        double blueWeight(size_t row, size_t col) const;
        double redSubWeight(Subgrid const& sg) const;
        double blueSubWeight(Subgrid const& sg) const;

        double yCoord(size_t row) const;
        double xCoord(size_t col) const;

        size_t size() const;
        Rectangle toRectangle(Subgrid const &sg) const;
    };

    inline Grid make_exact_grid(const wpoint_list_t& m_pts, const wpoint_list_t& b_pts) {
        return Grid(m_pts, b_pts);
    }

    inline Grid make_net_grid(const point_list_t& pts, const wpoint_list_t& m_pts, const wpoint_list_t& b_pts) {
        return Grid(pts, m_pts, b_pts);
    }


    Subgrid max_subgrid_convex(Grid const &grid, double eps, discrepancy_func_t const &f);
    Subgrid max_subgrid_linear(Grid const &grid, double a, double b);
    Subgrid max_subgrid(Grid const &grid, discrepancy_func_t const &func);

    //////////////////////////////////////////////////////////////
    /////Max Labeled rectangle code///////////////////////////////
    //////////////////////////////////////////////////////////////

    std::tuple<Rectangle, double> max_rect_labeled(size_t r, double max_w, lpoint_list_t const& m_points, lpoint_list_t const& b_points, const discrepancy_func_t& func);

    std::tuple<Rectangle, double> max_rect_labeled_scale(
            size_t r,
            double max_r,
            double alpha,
            const point_list_t &net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            const discrepancy_func_t &f);



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //FAST RECTANGLE SCANNING CODE//////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class EPoint {
        /*
         * We convert all of our points into Exact Points.
         */
    public:
        EPoint() : coord{0, 0}, weight(0.0) {}

        EPoint(size_t x, size_t y, double weight) : coord{x, y}, weight(weight) {}

        double get_weight() const { return weight; }
        void set_weight(double w) { weight = w; }
        size_t get_x() const { return coord[0]; }
        size_t get_y() const { return coord[1]; }

        bool operator<(const EPoint &b) const {
            //By default we sort by x value.
            return this->get_x() < b.get_x();
        }
        size_t operator()(const size_t& i) const {
            assert(i == 0 || i == 1);
            return coord[i];
        }

        size_t& operator()(const size_t& i) {
            assert(i == 0 || i == 1);
            return coord[i];
        }


    private:
        size_t coord[2];
        double weight;
    };

    using ept_t = EPoint;
    using epoint_list_t = std::vector<EPoint>;
    using epoint_it_t = std::vector<EPoint>::iterator;

    std::tuple<epoint_list_t, epoint_list_t> to_epoints(wpoint_list_t const& mpts, wpoint_list_t const& bpts);

    class ERectangle {
        size_t u_x, u_y, l_x, l_y;
    public:

        ERectangle(const ept_t& p1, const ept_t& p2, const ept_t& p3, const ept_t& p4) {
            u_x = std::max({p1(0), p2(0), p3(0), p4(0)});
            l_x = std::min({p1(0), p2(0), p3(0), p4(0)});
            u_y = std::max({p1(1), p2(1), p3(1), p4(1)});
            l_y = std::min({p1(1), p2(1), p3(1), p4(1)});
        }

        ERectangle() : u_x(0), u_y(0), l_x(0), l_y(0) {}

        ERectangle(size_t ux, size_t uy, size_t lx, size_t ly) : u_x(ux), u_y(uy), l_x(lx), l_y(ly) {}

        inline bool contains(const ept_t& p1) const {
            return (u_x > p1(0) && p1(0) >= l_x && u_y > p1(1) && p1(1) >= l_y);
        }

        std::string toString() const {
            std::stringstream ss;
            ss << "ERectangle(" << u_x << ", " << u_y << ", " << l_x << ", " << l_y << ")";
            return ss.str();
        }
        size_t lowX() const { return l_x; }
        size_t upX() const { return u_x; }
        size_t lowY() const { return l_y; }
        size_t upY() const { return u_y; }
    };


    template<typename R>
    double range_weight(const R& range, const epoint_list_t& pts) {
        double weight = 0.0;
        for (auto& pt: pts) {
            if (range.contains(pt)) {
                weight += pt.get_weight();
            }
        }
        return weight;
    }

    class Slab {

    public:
        std::vector<size_t> split_offsets;

        epoint_list_t m_merges;
        epoint_list_t b_merges;

        size_t top_y;
        size_t bottom_y;

        std::weak_ptr<Slab> parent;
        std::shared_ptr<Slab> up;
        std::shared_ptr<Slab> down;

        Slab(std::weak_ptr<Slab> p,
             epoint_list_t m_m,
             epoint_list_t b_m,
             size_t ty,
             size_t by);

        size_t get_mid() const {
            size_t midpoint;
            if (down != nullptr) midpoint = down->top_y;
            else if (up != nullptr) midpoint = up->bottom_y;
            else {
                midpoint = (top_y + bottom_y) / 2;
            }
            return midpoint;
        }

        bool has_mid() const {
            return down != nullptr || up != nullptr;
        }

        double measure_interval(size_t mxx, size_t mnx, double a, double b) const;
    };

    using slab_ptr = std::shared_ptr<Slab>;
    using wslab_ptr = std::weak_ptr<Slab>;
    using slab_it_t = std::vector<Slab>::const_iterator;

    class SlabTree {

    public:

        slab_ptr get_root() {
            return root;
        }

        slab_ptr get_containing(ERectangle const &rect) const;
        SlabTree(epoint_list_t ms, epoint_list_t bs, double max_w);
        SlabTree(std::vector<size_t> const &vert_decomp, epoint_list_t ms, epoint_list_t bs, bool compression, double max_w);

        //Useful for debuging the structure so you can define a non compressed version with fixed decomposition.
        void init(std::vector<size_t> const& vert_decomp, bool compression, double max_w);

        double measure_rect(ERectangle const &rect, double a, double b) const;
        std::tuple<ERectangle, double> max_rectangle(double m_a, double b_b);
    private:
        slab_ptr root;
        epoint_list_t mpts;
        epoint_list_t bpts;
    };


    std::tuple<Rectangle, double> max_rectangle(const wpoint_list_t& m_points, const wpoint_list_t& b_points, double eps, double a, double b);

}
#endif //PYSCAN_RECTANGLESCAN_HPP
