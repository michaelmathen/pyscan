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

        friend std::ostream& operator<<(std::ostream& os, EPoint const& el) {
            os << "ept({" << el.coord[0] << ", "  << el.coord[1] << "}, "<<  el.weight << ")";
            return os;
        }

    private:
        size_t coord[2];
        double weight;
    };

    using ept_t = EPoint;
    using epoint_list_t = std::vector<EPoint>;
    using epoint_it_t = std::vector<EPoint>::iterator;

    std::tuple<epoint_list_t, epoint_list_t, std::unordered_map<size_t, double>, std::unordered_map<size_t, double> >
            to_epoints(wpoint_list_t const& mpts, wpoint_list_t const& bpts);

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

    inline double computeTotal(const std::vector<EPoint>& pts) {
        double res = 0.0;
        for (auto& x: pts) res += x.get_weight();
        return res;
    }

    template <typename R>
    double evaluate_range(
            const R& range,
            const std::vector<EPoint>& red,
            const std::vector<EPoint>& blue,
            const discrepancy_func_t& f) {

        return f(range_weight(range, red), computeTotal(red),
                 range_weight(range, blue), computeTotal(blue));
    }

    template <typename R>
    std::tuple<R, double> max_range4(
            const std::vector<EPoint>& vertical_point_net,
            const std::vector<EPoint>& horz_point_net,
            const std::vector<EPoint>& red,
            const std::vector<EPoint>& blue,
            const discrepancy_func_t& f) {

        double max_stat = 0.0;
        R cur_max;
        for (size_t i = 0; i < vertical_point_net.size() - 1; ++i) {
            for (size_t j = i + 1; j < vertical_point_net.size(); ++j) {
                for (size_t k = 0; k < horz_point_net.size() - 1; ++k) {
                    for (size_t l = k + 1; l < horz_point_net.size(); ++l) {
                        auto y_max = std::max(vertical_point_net[i].get_y(), vertical_point_net[j].get_y());
                        auto y_min = std::min(vertical_point_net[i].get_y(), vertical_point_net[j].get_y());
                        auto x_max = std::max(horz_point_net[i].get_x(), horz_point_net[j].get_x());
                        auto x_min = std::min(horz_point_net[i].get_x(), horz_point_net[j].get_x());

                        R now(x_max, y_max, x_min, y_min);
                        double cur_stat = evaluate_range(now, red, blue, f);
                        if (cur_stat > max_stat) {
                            cur_max = now;
                            max_stat = cur_stat;
                        }
                    }
                }
            }
        }

        return std::make_tuple(cur_max, max_stat);
    }

    class Interval {
        size_t left;
        size_t right;
        double value;
    public:
        Interval(size_t l, size_t r, double v) : left(l), right(r), value(v) {}

        size_t get_r() const { return right; }
        size_t get_l() const { return left; }
        double_t get_v() const {return value; }

        friend std::ostream& operator<<(std::ostream& os, Interval const& el) {
            os << "Interval(" << el.left << ", " <<  el.right << ", " <<  el.value << ")";
            return os;
        }
    };


    class MaxIntervalAlt {
        Interval left_max;
        Interval right_max;
        Interval center_max;
    public:
        MaxIntervalAlt(size_t lpt, size_t rpt);
        MaxIntervalAlt(size_t val, double weight);
        MaxIntervalAlt(Interval l, Interval c, Interval r) : left_max(l), right_max(r), center_max(c) {}
        MaxIntervalAlt &operator+=(const MaxIntervalAlt& op);
//         MaxIntervalAlt &operator+=(const MaxIntervalAlt& op) {
//             right_max = op.right_max;
//             return *this;
//         }

        friend std::ostream& operator<<(std::ostream& os, MaxIntervalAlt const& el) {
            os << "MaxIntervalAlt(" << el.left_max << ", "  << el.center_max << ", "<< el.right_max << "]";
            return os;
        }

        //void update_left_weight(double weight);

        size_t left() const;
        size_t right() const;

        Interval get_max() const;
    };

    class Slab {

    public:
        //This is the union of m_merges and b_merges and all child m_merges and b_merges.
        std::vector<size_t> global_split_offset;

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

        friend std::ostream& operator<<(std::ostream& os, std::shared_ptr<Slab> el) {
            if (el == nullptr) {
                os << "[]";
            } else {
                os << el->global_split_offset;
            }
            return os;
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
        SlabTree(slab_ptr child, double tm, double tb) : root(child), total_m(tm), total_b(tb) {}

        SlabTree(epoint_list_t ms, epoint_list_t bs, double max_w);
        SlabTree(std::vector<size_t> const &vert_decomp, epoint_list_t ms, epoint_list_t bs, bool compression, double max_w);

        //Useful for debuging the structure so you can define a non compressed version with fixed decomposition.
        void init(epoint_list_t mpts, epoint_list_t bpts, std::vector<size_t> const &vert_decomp, bool compression, double max_w);

        double measure_rect(ERectangle const &rect, double a, double b) const;
        std::tuple<ERectangle, double> max_rectangle(double m_a, double b_b);

        friend std::ostream& operator<<(std::ostream& os, SlabTree const& el) {
            std::queue<std::tuple<size_t, decltype(root)>> curr_queue;
            curr_queue.emplace(0, el.root);
            os << "TREE" << std::endl;
            os << el.root->global_split_offset << std::endl;
            os << "0 = ";
            size_t curr_level = 0;
            while (!curr_queue.empty()) {
                auto [l, el] = curr_queue.front();
                if (l != curr_level) {
                    os << std::endl;
                    os << l << " = ";
                    curr_level = l;
                }
                os << "[" << el->top_y << ", " << el->bottom_y << "]";
                curr_queue.pop();
                if (el->up != nullptr) {
                    curr_queue.emplace(l + 1, el->up);
                }
                if (el->down != nullptr) {
                    curr_queue.emplace(l + 1, el->down);
                }
            }
            return os;
        }

        SlabTree get_upper_tree() const {
            return SlabTree(root->up, total_m, total_b);
        }

        SlabTree get_lower_tree() const {
            return SlabTree(root->down, total_m, total_b);
        }
    private:
        slab_ptr root;
        double total_m;
        double total_b;
    };

    std::vector<MaxIntervalAlt> insert_updates(std::vector<MaxIntervalAlt> const& max_intervals,
                                               epoint_list_t const& updates, double scale);

    std::vector<MaxIntervalAlt> reduce_merges(std::vector<MaxIntervalAlt> const& max_intervals,
                                              std::vector<size_t> const& curr_splits);

    std::tuple<Rectangle, double> max_rectangle(const wpoint_list_t& m_points, const wpoint_list_t& b_points, double eps, double a, double b);




    class LeafNode {
    public:

        LeafNode(size_t pos, double weight) : max_interval(pos, weight) {}

        virtual size_t get_left_bound() const {
            return max_interval.left();
        }

        virtual size_t get_right_bound() const {
            return max_interval.right();
        }

        virtual std::shared_ptr<LeafNode> get_right() const {
             return nullptr;
        }

        virtual  std::shared_ptr<LeafNode> get_left() const {
            return nullptr;
        }

        virtual bool contains(size_t x) {
             return get_left_bound() <= x && x < get_right_bound();
         }


        virtual void rebuild() {
        }

        virtual size_t get_mid() const {
            return  get_left_bound();
        }

        virtual bool is_leaf() const { return true; }

        virtual void set_parent(std::shared_ptr<LeafNode> lp) {
            parent = lp;
        }

        virtual bool is_left(std::shared_ptr<LeafNode> ptr) const { return false; }

    protected:
        friend class IntervalTreap;
        friend class IntervalNode;
        MaxIntervalAlt max_interval;
        std::weak_ptr<LeafNode> parent;
    };

    class IntervalNode : public LeafNode {

    public:

        IntervalNode(std::shared_ptr<LeafNode> l1, std::shared_ptr<LeafNode> l2, double priority) : LeafNode(*l1),
            r_child(l2), l_child(l1), priority(priority) {
            max_interval += l2->max_interval;
        }

        bool is_leaf() const override { return false; }

        size_t get_mid() const override {
            return r_child->get_left_bound();
        }

        void rebuild() override {
            /*
             * Builds this node from its children by merging the left child and right child max intervals.
             */
            max_interval = l_child->max_interval;
            max_interval += r_child->max_interval;
        }

        friend std::shared_ptr<LeafNode> l_rotation(std::shared_ptr<IntervalNode> p) {
            /*
             *     p               r
             *  l    r     =>   p    rr
             *      lr rr     l  lr
             *
             */
            auto r = p->r_child;
            auto l = p->l_child;
            auto lr = p->r_child->get_left();

            assert(!r->is_leaf());
            ((IntervalNode*)r.get())->l_child = p;
            p->r_child = lr;

            r->set_parent(p->parent.lock());
            p->set_parent(r);
            lr->set_parent(p);

            p->rebuild();
            l->rebuild();
            return r;
        }

        friend std::shared_ptr<LeafNode> r_rotation(std::shared_ptr<IntervalNode> p) {
            /*
             *     p               l
             *  l     r     =>   ll    p
             * ll rl                rl   r
             */
            auto r = p->r_child;
            auto l = p->l_child;
            auto rl = p->l_child->get_right();

            assert(!l->is_leaf());
            ((IntervalNode*)l.get())->r_child = p;
            p->l_child = rl;

            l->set_parent(p->parent.lock());
            p->set_parent(l);
            rl->set_parent(p);

            p->rebuild();
            l->rebuild();
            return l;
        }

        bool is_left(std::shared_ptr<LeafNode> ptr) const override { return ptr == l_child; }

    protected:
        friend class IntervalTreap;
        std::shared_ptr<LeafNode> r_child;
        std::shared_ptr<LeafNode> l_child;
        double priority;
    };

    class IntervalTreap {
    public:

        IntervalTreap() : root(nullptr), distribution(0.0, 1.0), gen(std::random_device()()) {}

        std::shared_ptr<LeafNode> upper_bound(size_t val) const {
            /*
             * Finds the lowest node that we are contained inside.
             */
            auto curr_root = root;
            auto p = curr_root;
            while (curr_root != nullptr) {
                p = curr_root;
                if (curr_root->get_mid() > val) {
                    curr_root = curr_root->get_left();
                } else {
                    curr_root = curr_root->get_right();
                }
            }
            return p;
        }

        void insert(size_t val, double weight) {
            auto c_right = upper_bound(val);
            std::shared_ptr<LeafNode> c_left = std::make_shared<LeafNode>(val, weight);

            if (c_right == nullptr) {
                root = c_left;
            } else {
                auto p = c_right->parent.lock();
                std::shared_ptr<LeafNode>* child = &root;
                if (p != nullptr) {
                    if (p->is_left(c_right)) {
                        child = &(std::dynamic_pointer_cast<IntervalNode>(p)->l_child);
                    } else {
                        child = &(std::dynamic_pointer_cast<IntervalNode>(p)->r_child);
                    }
                }
                if (c_right->get_mid() <= val) {
                    std::swap(c_left, c_right);
                }
                *child = std::make_shared<LeafNode>(new IntervalNode(c_left, c_right, distribution(gen)));
                (*child)->set_parent(p);

                // Binary tree now.

                //Now have to bubble up the change to the root.
                while ((*child)->parent.lock() != nullptr) {
                    auto p = std::dynamic_pointer_cast<IntervalNode>((*child)->parent.lock());
                    auto curr = std::dynamic_pointer_cast<IntervalNode>(*child);

                    if (p->parent.lock() == nullptr) {
                        child = &root;
                    } else {
                        if (p->parent.lock()->is_left(p)) {
                            child = &(std::dynamic_pointer_cast<IntervalNode>(p)->l_child);
                        } else {
                            child = &(std::dynamic_pointer_cast<IntervalNode>(p)->r_child);
                        }
                    }

                    if (p->priority > curr->priority) {
                        if (p->is_left(curr)) {
                            *child = r_rotation(p);
                        } else {
                            *child = l_rotation(p);
                        }
                    }

                }

            }
        }

    private:
        std::shared_ptr<LeafNode> root;
        std::uniform_real_distribution<double> distribution;
        std::minstd_rand gen;

    };


    std::tuple<ERectangle, double> max_rectangle_heap(epoint_list_t const& mpts, epoint_list_t const& bpts, double x, double y) {

    }
}
#endif //PYSCAN_RECTANGLESCAN_HPP
