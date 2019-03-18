//
// Created by mmath on 3/3/19.
//

#include "Utilities.hpp"
#include "IntervalScan.hpp"

namespace pyscan {

    MaxIntervalAlt::MaxIntervalAlt(size_t val, double weight) : Interval(val, val, weight), left_max(val, val, weight), right_max(val, val, weight),
                                                                center_max(val, val, weight) {}

    MaxIntervalAlt::MaxIntervalAlt(size_t lpt, size_t rpt) : Interval(lpt, rpt, 0.0), left_max(lpt, rpt, 0.0), right_max(lpt, rpt, 0.0),
                                                             center_max(lpt, rpt, 0.0), full_weight(0.0) {}



    Interval arg_max(std::initializer_list<Interval> args) {
        Interval max_el = *args.begin();
        for (auto el = args.begin() + 1; args.end() != el; ++el) {
            if (el->get_v() > max_el.get_v()) {
                max_el = *el;
            }
        }
        return max_el;
    }

    MaxIntervalAlt &MaxIntervalAlt::operator+=(const MaxIntervalAlt &op) {

        center_max = arg_max({center_max, op.center_max, right_max + op.left_max});

        right_max = arg_max({op.right_max, right_max + op});
        left_max = arg_max({left_max, *this + op.left_max});

        this->right = op.get_r();
        this->value = op.get_v() + get_v();
        return *this;
    }

    Interval MaxIntervalAlt::get_max() const {
        return center_max;
    }


    Interval max_interval(std::vector<size_t> indices, std::vector<double> weights) {
        /*
         * Finds the max interval over this set of weighted indices.
         */
        if (indices.empty()) {
            return Interval(0, 0, 0.0);
        }
        auto perm = util::sort_permutation(indices, std::less<>());
        util::apply_permutation_in_place(indices, perm);
        util::apply_permutation_in_place(weights, perm);

        size_t begin = 0;
        Interval max_interval(0, 0, weights[0]);
        double max_ending_here = weights[0];
        for (size_t i = 1; i < indices.size(); ++i) {
            if (max_ending_here + weights[i] <= weights[i]) {
                max_ending_here = weights[i];
                begin = indices[i];
            } else {
                max_ending_here += weights[i];
            }

            if (max_interval.get_v() <= max_ending_here) {
                max_interval = Interval(begin, indices[i], max_ending_here);
            }
        }
        return max_interval;
    }

    Interval max_interval_slow(std::vector<size_t> indices, std::vector<double> weights) {
        /*
         * Finds the max interval over this set of weighted indices.
         */
        if (indices.empty()) {
            return Interval(0, 0, 0.0);
        }
        auto perm = util::sort_permutation(indices, std::less<>());
        util::apply_permutation_in_place(indices, perm);
        util::apply_permutation_in_place(weights, perm);
        Interval max_interval(0, 0, weights[0]);

        for (size_t i = 0; i < indices.size(); i++) {
            double w = 0;
            for (size_t j = i; j < indices.size(); j++) {
                w += weights[j];
                if (w > max_interval.get_v()) {
                    max_interval = Interval(indices[i], indices[j], w);
                }
            }
        }
        return max_interval;
    }
}