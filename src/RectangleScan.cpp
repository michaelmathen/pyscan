//
// Created by mmath on 5/28/17.
//
#include <functional>
#include <tuple>

#include "RectangleScan.hpp"

namespace pyscan {


     void mergeWeights(MaximumIntervals &interval, slab_ptr left, slab_ptr right) {
         interval = interval.mergeZeros(left->non_zeros, right->non_zeros);
     }

    void updateWeights(MaximumIntervals &interval, slab_ptr &slab, double a, double b) {
        interval.updateWeights(slab->red_weights, slab->indices, a);
        interval.updateWeights(slab->blue_weights, slab->indices, b);
        interval.updateBounds(slab->left_col, slab->right_col);
    }


    Subgrid maxSubgridLinearLeft(slab_ptr slab, MaximumIntervals const& maxInt, double a, double b) {
        if (slab->hasChildren()) {
            MaximumIntervals case1 = maxInt;
            updateWeights(case1, slab->right, a, b);
            case1 = case1.mergeZeros(slab->left->non_zeros);
            auto case2 = maxInt.mergeZeros(slab->right->non_zeros);
            auto s1 = maxSubgridLinearLeft(slab->left, case1, a, b);
            auto s2 = maxSubgridLinearLeft(slab->right, case2, a, b);
            return s1.fValue() > s2.fValue() ? s1 : s2;
        } else {
            MaximumIntervals case1 = maxInt;
            updateWeights(case1, slab, a, b);
            return case1.getMax();
        }
    }

    Subgrid maxSubgridLinearRight(slab_ptr slab, MaximumIntervals const& maxInt, double a, double b) {
        if (slab->hasChildren()) {
            MaximumIntervals case1 = maxInt;
            updateWeights(case1, slab->left, a, b);
            case1 = case1.mergeZeros(slab->right->non_zeros);
            auto case2 = maxInt.mergeZeros(slab->left->non_zeros);
            auto s1 = maxSubgridLinearRight(slab->right, case1, a, b);
            auto s2 = maxSubgridLinearRight(slab->left, case2, a, b);
            return s1.fValue() > s2.fValue() ? s1 : s2;
        } else {
            MaximumIntervals case1 = maxInt;
            updateWeights(case1, slab, a, b);
            return case1.getMax();
        }
    }

    Subgrid maxSubgridSpan(slab_ptr left_side, slab_ptr right_side, MaximumIntervals const &maxInt,
                             double a, double b) {

        using SubProblem = std::tuple<slab_ptr, slab_ptr, MaximumIntervals>;
        std::deque<SubProblem> sub_problem_stack;
        sub_problem_stack.push_back(make_tuple(left_side, right_side, maxInt));
        Subgrid max_subgrid = Subgrid(-1, -1, -1, -1, -std::numeric_limits<double>::infinity());

        while (sub_problem_stack.size() > 0) {
           SubProblem sub_problem = sub_problem_stack.back();
           sub_problem_stack.pop_back();
           auto& left = std::get<0>(sub_problem);
           auto& right = std::get<1>(sub_problem);
           auto ll = left->left;
           auto lr = left->right;
           auto rl = right->left;
           auto rr = right->right;
           Subgrid new_max = Subgrid(-1, -1, -1, -1, -std::numeric_limits<double>::infinity());
           if (!left->hasChildren() && !right->hasChildren()) {
              //std::cout << "no Children" << std::endl;
              MaximumIntervals case1 = maxInt;
              updateWeights(case1, left, a, b);
              updateWeights(case1, right, a, b);
              new_max = case1.getMax();
           } else if (!left->hasChildren()) {
               //std::cout << "left case" << std::endl;
               MaximumIntervals case1 = maxInt;
               updateWeights(case1, left, a, b);
               new_max = maxSubgridLinearRight(right, case1, a, b);
           } else if (!right->hasChildren()) {
               //std::cout << "right case" << std::endl;
               MaximumIntervals case1 = maxInt;
               updateWeights(case1, right, a, b);
               new_max = maxSubgridLinearLeft(left, case1, a, b);
           } else {
                SubProblem  sub1 = SubProblem(ll, rr, maxInt);
                auto& case1 = std::get<2>(sub1);
                updateWeights(case1, lr, a, b);
                updateWeights(case1, rl, a, b);
                mergeWeights(case1, ll, rr);

                SubProblem  sub2 = SubProblem(ll, rl, maxInt);
                auto& case2 = std::get<2>(sub2);
                updateWeights(case2, lr, a, b);
                mergeWeights(case2, ll, rl);

                SubProblem  sub3 = SubProblem(lr, rr, maxInt);
                auto& case3 = std::get<2>(sub3);
                updateWeights(case3, rl, a, b);
                mergeWeights(case3, lr, rr);

                SubProblem sub4 = SubProblem(lr, rl, maxInt);
                auto& case4 = std::get<2>(sub4);
                mergeWeights(case4, lr, rl);
            }
            if (new_max.fValue() > max_subgrid.fValue()) {
                  max_subgrid = new_max;
            }
        }
        return max_subgrid;
    }


    Subgrid maxSubgridSpanRec(slab_ptr left, slab_ptr right, MaximumIntervals const &maxInt,
                             double a, double b) {
        //std::cout << maxInt << maxInt.getLeft() << "  " << maxInt.getRight() <<  std::endl;
        //std::cout << maxInt << std::endl;
        //std::cout << maxInt.getIntervalNum() << " " << maxInt.getWeightNum() << std::endl;
        //std::cout << maxInt << std::endl;


        auto ll = left->left;
        auto lr = left->right;
        auto rl = right->left;
        auto rr = right->right;

        if (!left->hasChildren() && !right->hasChildren()) {
            //std::cout << "no Children" << std::endl;
            MaximumIntervals case1 = maxInt;
            updateWeights(case1, left, a, b);
            updateWeights(case1, right, a, b);
            return case1.getMax();
        }

        if (!left->hasChildren()) {
            //std::cout << "left case" << std::endl;
            MaximumIntervals case1 = maxInt;
            updateWeights(case1, left, a, b);
            return maxSubgridLinearRight(right, case1, a, b);
        }

        if (!right->hasChildren()) {
            //std::cout << "right case" << std::endl;
            MaximumIntervals case1 = maxInt;
            updateWeights(case1, right, a, b);
            return maxSubgridLinearLeft(left, case1, a, b);
        }

        //std::cout << "Standard case" << std::endl;
        MaximumIntervals case1 = maxInt;
        MaximumIntervals case2 = maxInt;
        MaximumIntervals case3 = maxInt;
        MaximumIntervals case4 = maxInt;

        updateWeights(case1, lr, a, b);
        updateWeights(case1, rl, a, b);
        mergeWeights(case1, ll, rr);
        auto s1 = maxSubgridSpan(ll, rr, case1, a, b);
        updateWeights(case2, lr, a, b);
        mergeWeights(case2, ll, rl);
        auto s2 = maxSubgridSpan(ll, rl, case2, a, b);
        updateWeights(case3, rl, a, b);
        mergeWeights(case3, lr, rr);
        auto s3 = maxSubgridSpan(lr, rr, case3, a, b);
        mergeWeights(case4, lr, rl);
        auto s4 = maxSubgridSpan(lr, rl, case4, a, b);
        return std::max({s1, s2, s3, s4}, [](Subgrid const &e1, Subgrid const &e2) {
            return e1.fValue() < e2.fValue();
        });
    }

    Subgrid maxSubgridLinear(slab_ptr slab, MaximumIntervals const &maxInt, double a, double b) {
        if (slab == nullptr) {
            //std::cout << maxInt << maxInt.getLeft() << maxInt.getRight() <<  std::endl;
            return maxInt.getMax();
        } else {
            Subgrid max = Subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
            std::deque<std::tuple<slab_ptr, MaximumIntervals>> slabs;
            slabs.push_back(std::make_tuple(slab, maxInt));
            while (slabs.size() > 0) {
                auto slab_mx = slabs.back();
                auto new_slab = std::get<0>(slab_mx);
                auto mx = std::get<1>(slab_mx);
                slabs.pop_back();
                auto new_maxInt = mx.mergeZeros(new_slab->non_zeros);
                if (new_slab->left == nullptr && new_slab->right == nullptr) {
                    updateWeights(new_maxInt, new_slab, a, b);
                    auto new_max = new_maxInt.getMax();
                    if (new_max.fValue() > max.fValue()) {
                        max = new_max;
                    }
                    continue;
                } else {
                    new_maxInt.setBounds(new_slab->getMid() + 1, new_slab->getMid());
                    auto new_max = maxSubgridSpan(new_slab->left, new_slab->right, new_maxInt, a, b);
                    if (new_max.fValue() > max.fValue()) {
                        max = new_max;
                    }
                    slabs.push_back(make_tuple(new_slab->left, new_maxInt));
                    slabs.push_back(make_tuple(new_slab->right, new_maxInt));
                }
            }
            return max;
        }
    }
}

