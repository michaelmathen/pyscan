//
// Created by mmath on 5/28/17.
//
#include <functional>
#include <tuple>
#include <unordered_set>

#include "RectangleScan.hpp"

namespace pyscan {

    /*
     * Simple 1/eps^4 algorithm described in the paper that just computes every subgrid.
     * This will work on a nonlinear function.
     */
    template<typename Weight, typename F>
    Subgrid maxSubgridNonLinear(Grid<Weight> const &grid, F func, double rho) {
        std::vector<Weight> red_count(grid.size(), 0);
        std::vector<Weight> blue_count(grid.size(), 0);
        double t_red = static_cast<double>(grid.totalRedWeight());
        double t_blue = static_cast<double>(grid.totalBlueWeight());
        Subgrid max = Subgrid(-1, -1, -1, -1, -std::numeric_limits<double>::infinity());
        for (size_t i = 0; i < grid.size() - 1; i++) {
            red_count.assign(grid.size(), 0);
            blue_count.assign(grid.size(), 0);
            for (size_t j = i; j < grid.size(); j++) {

                for (size_t k = 0; k < grid.size(); k++) {
                    blue_count[k] += grid.blueCount(j, k);
                    red_count[k] += grid.redCount(j, k);
                }

                for (size_t k = 0; k < grid.size(); k++) {
                    Weight red_l_count = 0;
                    Weight blue_l_count = 0;
                    for (size_t l = k; l < grid.size(); l++) {
                        red_l_count += red_count[l];
                        blue_l_count += blue_count[l];
                        double maxf = func(red_l_count / t_red, blue_l_count / t_blue, rho);
                        //std::cout << red_l_count / t_red << " " << blue_l_count / t_blue << " " << maxf << std::endl;
                        if (maxf > max.fValue()) {
                            max = Subgrid(l, j, k, i, maxf);
                        }
                    }

                }
            }
        }
        return max;
    }


    Subgrid maxSubgridKullSlow(Grid<double> const &grid, double rho) {
        return maxSubgridNonLinear(grid, kulldorff, rho);
    }

    Subgrid maxSubgridLinearSlow(Grid<int> const& grid, double a, double b) {
        return maxSubgridNonLinear(grid, [&](double red, double blue, double rho) {
            return a * red + b * blue;
        }, 0);
    }

    Subgrid maxSubgridKullSlow(Grid<int> const &grid, double rho) {
        return maxSubgridNonLinear(grid, kulldorff, rho);
    }

    Subgrid maxSubgridLinearSlow(Grid<double> const& grid, double a, double b) {
        return maxSubgridNonLinear(grid, [&](double red, double blue, double rho) {
            return a * red + b * blue;
        }, 0);
    }

    template<typename Weight>
    Subgrid maxSubgridLinearSimpleInt(Grid<Weight> const &grid, double a, double b) {
        std::vector<double> weight(grid.size(), 0);
        Subgrid max = Subgrid(-1, -1, -1, -1, -std::numeric_limits<double>::infinity());
        for (size_t i = 0; i < grid.size(); i++) {
            weight.assign(grid.size(), 0);
            for (size_t j = i; j < grid.size(); j++) {
                for (size_t k = 0; k < grid.size(); k++) {
                    weight[k] += b * grid.blueWeight(j, k);
                    weight[k] += a * grid.redWeight(j, k);
                }
                // This is just computing the max interval over the weight vector.
                // Basically we scan until the weight drops to less than 0 at which point
                // we can just restart the interval here and it will always be larger.
                double curr_W = 0;
                size_t start_ix = 0;
                for (size_t l = 0; l < grid.size(); l++) {
                    curr_W += weight[l];
                    if (curr_W <= 0) {
                        curr_W = 0;
                        start_ix = l + 1;
                    }
                    if (curr_W > max.fValue()) {
                        max = Subgrid(l, j, start_ix, i, curr_W);
                    }
                }
            }
        }
        return max;
    }

    Subgrid maxSubgridLinearSimple(Grid<double> const& grid, double a, double b) {
        return maxSubgridLinearSimpleInt(grid, a, b);
    }

    Subgrid maxSubgridLinearSimple(Grid<int> const& grid, double a, double b) {
        return maxSubgridLinearSimpleInt(grid, a, b);
    }
    //template <> Subgrid MaxSubgridLinearSimple<int>(Grid<int> const& grid, double a, double b);

    template <typename W>
    double getMeasured(Point<W, 2> const& pt) {
        return pt.getRedWeight();
    }

    template <typename W>
    double getBaseline(Point<W, 2> const& pt) {
        return pt.getBlueWeight();
    }

    template< typename T, typename F>
    double computeLabelTotal(T begin, T end, F func) {
        std::unordered_set<size_t> label_set;
        double total = 0;
        for (; begin != end; ++begin) {
            if (label_set.end() == label_set.find(begin->getLabel())) {
                total += func(*begin);
                label_set.insert(begin->getLabel());
            }
        }
        return total;
    }


    template <typename W, typename F>
    Rectangle maxLabeledRect(std::vector<LPoint<W>> const& net,
                             std::vector<LPoint<W>> const& m_points,
                             std::vector<LPoint<W>> const& b_points,
                             F func) {

        double m_Total = computeLabelTotal(m_points.begin(), m_points.end(), getMeasured<W>);
        double b_Total = computeLabelTotal(b_points.begin(), b_points.end(), getBaseline<W>);
        Rectangle maxRect(0, 0, 0, 0, 0.0);
        //Create a top and bottom line
        for (auto onb1 = net.begin(); onb1 != (net.end() - 2); onb1++) {
            for (auto onb2 = onb1 + 1; onb2 != (net.end() - 1); onb2++) {
                auto nb1 = get<1>(*onb1) >= get<1>(*onb2) ? onb1 : onb2;
                auto nb2 = get<1>(*onb1) < get<1>(*onb2) ? onb1 : onb2;

                std::vector<LPoint<W>> divisions;
                std::vector<LPoint<W>> m_slab;
                std::vector<LPoint<W>> b_slab;
                auto between_f = [&](LPoint<W> const& pt){
                    return get<1>(pt) <= get<1>(*nb1) && get<1>(*nb2) <= get<1>(pt);
                };
                auto pointXComp = [](LPoint<W> const& p1, LPoint<W> const& p2){
                    return get<0>(p1) < get<0>(p2);
                };
                std::remove_copy_if(net.begin(), net.end(), std::back_inserter(divisions), between_f);
                //Filter the m_slab and b_slab points
                std::remove_copy_if(m_points.begin(), m_points.end(), std::back_inserter(m_slab), between_f);
                std::remove_copy_if(b_points.begin(), b_points.end(), std::back_inserter(b_slab), between_f);
                std::sort(divisions.begin(), divisions.end(), pointXComp);

                //Now position the points into each division
                using part_t = std::vector<std::vector<LPoint<W>>>;
                part_t m_parts(divisions.size() + 1,  std::vector<LPoint<W>>());
                part_t b_parts(divisions.size() + 1, std::vector<LPoint<W>>());

                auto place_pts = [&] (part_t & parts, decltype(m_points) all_pts) {
                    for (auto& pt : all_pts) {
                        auto loc_it = std::lower_bound(divisions.begin(), divisions.end(), pt, pointXComp);
                        parts[loc_it - divisions.begin()].emplace_back(pt);
                    }
                };
                place_pts(m_parts, m_slab);
                place_pts(b_parts, b_slab);

                for (int i = 0; i < divisions.size() - 1; i++) {
                    double m_count = 0;
                    double b_count = 0;
                    std::unordered_set<size_t> active_m_labels;
                    std::unordered_set<size_t> active_b_labels;
                    for (int j = i + 1; j < divisions.size(); j++) {
                        for (auto m_b = m_parts[j].begin(); m_b != m_parts[j].end(); ++m_b) {
                            if (active_m_labels.find(m_b->getLabel()) != active_m_labels.end()) {
                                m_count += m_b->getRedWeight();
                                active_m_labels.insert(m_b->getLabel());
                            }
                        }
                        for (auto b_b = b_parts[j].begin(); b_b != b_parts[j].end(); ++b_b) {
                            if (active_m_labels.find(b_b->getLabel()) != active_b_labels.end()) {
                                b_count += b_b->getRedWeight();
                                active_b_labels.insert(b_b->getLabel());
                            }
                        }
                        double new_val = func(m_count / m_Total, b_count / b_Total);
                        if (new_val > maxRect.fValue()) {
                            maxRect = Rectangle(get<0>(divisions[i]), get<1>(*nb1), get<0>(divisions[j]), get<1>(*nb2), new_val);
                        }
                    }
                }
            }
        }
        return maxRect;
    }

    Rectangle maxLabeledRectStat(std::vector<LPoint<int>> const& net,
                                 std::vector<LPoint<int>> const& m_points,
                                 std::vector<LPoint<int>> const& b_points, double rho) {
        return maxLabeledRect(net, m_points, b_points, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }

    Rectangle maxLabeledRectStat(std::vector<LPoint<double>> const& net,
                                 std::vector<LPoint<double>> const& m_points,
                                 std::vector<LPoint<double>> const& b_points, double rho) {
        return maxLabeledRect(net, m_points, b_points, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }

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
