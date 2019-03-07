//
// Created by mmath on 3/3/19.
//

#ifndef PYSCAN_INTERVALSCAN_HPP
#define PYSCAN_INTERVALSCAN_HPP

#include <cassert>
#include <memory>
#include <vector>

namespace pyscan {
    class Interval {
    protected:
        size_t left;
        size_t right;
        double value;
    public:
        Interval(size_t l, size_t r, double v) : left(l), right(r), value(v) {}

        virtual size_t get_r() const { return right; }
        virtual size_t get_l() const { return left; }
        virtual double_t get_v() const {return value; }

        friend Interval operator+(Interval const& e1, Interval const& e2) {
            return Interval(e1.get_l(), e2.get_r(), e1.get_v() + e2.get_v());
        }

        friend std::ostream& operator<<(std::ostream& os, Interval const& el) {
            os << "Interval(" << el.left << ", " <<  el.right << ", " <<  el.value << ")";
            return os;
        }
    };


    class MaxIntervalAlt : public Interval {
        Interval left_max;
        Interval right_max;
        Interval center_max;
        double full_weight;
    public:
        MaxIntervalAlt(size_t lpt, size_t rpt);
        MaxIntervalAlt(size_t val, double weight);
        //MaxIntervalAlt(Interval l, Interval c, Interval r) : left_max(l), right_max(r), center_max(c) {}
        MaxIntervalAlt &operator+=(const MaxIntervalAlt& op);

        friend std::ostream& operator<<(std::ostream& os, MaxIntervalAlt const& el) {
            os << "MaxIntervalAlt(" << el.left_max << ", "  << el.center_max << ", "<< el.right_max << "]";
            return os;
        }

        //void update_left_weight(double weight);


        Interval get_max() const;
    };


    class LeafNode {
    public:

        LeafNode(size_t pos, double weight) : max_interval(pos, weight) {}

        virtual size_t get_left_bound() const {
            return max_interval.get_l();
        }

        virtual size_t get_right_bound() const {
            return max_interval.get_r();
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

        virtual bool is_left(std::shared_ptr<LeafNode> ptr) const {
            (void)ptr;
            return false;
        }

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
                *child = std::dynamic_pointer_cast<LeafNode>(std::make_shared<IntervalNode>(c_left, c_right, distribution(gen)));
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


    Interval max_interval(std::vector<size_t> indices, std::vector<double> weights);

    Interval max_interval_slow(std::vector<size_t> indices, std::vector<double> weights);
}
#endif //PYSCAN_INTERVALSCAN_HPP
