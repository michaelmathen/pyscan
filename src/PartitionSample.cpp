//
// Created by mmath on 12/9/18.
//

#include <queue>
#include <random>
#include "Sampling.hpp"
#include "HalfSpaceScan.hpp"
#include "PartitionSample.hpp"


namespace pyscan {

    using split_node_t = std::tuple<wpoint_list_t, wpoint_list_t>;
    using split_queue_t = std::vector<split_node_t>;

    using lsplit_node_t = std::tuple<lpoint_list_t, lpoint_list_t>;
    using lsplit_queue_t = std::vector<split_node_t>;


    split_node_t split_directional(wpoint_list_t pts, double x, double y) {
        wpoint_list_t upper;
        wpoint_list_t lower;
        std::sort(pts.begin(), pts.end(), [x, y] (const pt2_t& p1, const pt2_t& p2) {
            return p1(0) * x  + p1(1) * y  < p2(0) * x  + p2(1) * y;
        });
        lower.insert(lower.end(), pts.begin(), pts.begin() + pts.size() / 2);
        upper.insert(upper.end(), pts.begin() + pts.size() / 2, pts.end());
        return std::make_tuple(upper, lower);
    }

    wpoint_list_t ham_tree_sample(const wpoint_list_t &pts, size_t s_size) {

        if (pts.size() <= s_size) {
            return pts;
        }
        auto cmp = [](const split_node_t& lhs, const split_node_t& rhs) {
            auto& [lu, ll] = lhs;
            auto& [ru, rl] = rhs;
            return lu.size() + ll.size() < ru.size() + rl.size();
        };

        std::random_device rd;
        std::minstd_rand gen(rd());
        split_queue_t split_queue = {split_directional(pts, 0, 1)};

        while (!split_queue.empty() && split_queue.size() < s_size / 2) {
            pop_heap(split_queue.begin(), split_queue.end(), cmp);
            auto [upper, lower] = split_queue[split_queue.size() - 1];
            split_queue.pop_back();
            //Take a net sample.

            //Add the wf back.
            auto net_sample1 = random_sample_wor(upper, gen, HAM_NET_SIZE / 2);
            auto net_sample2 = random_sample_wor(lower, gen, HAM_NET_SIZE / 2);
            point_list_t net_pts;
            for (auto& p : net_sample1) {
                net_pts.emplace_back(p);
            }
            for (auto& p : net_sample2) {
                net_pts.emplace_back(p);
            }

            if (upper.empty()) {
                split_queue.push_back(split_directional(upper, 0, 1));
                std::push_heap(split_queue.begin(), split_queue.end(), cmp);
            } else if (lower.empty()) {
                split_queue.push_back(split_directional(lower, 0, 1));
                std::push_heap(split_queue.begin(), split_queue.end(), cmp);
            } else {

                auto[plane, hTotal] = max_halfplane(net_pts, upper, lower,
                                                    [](double m, double m_total, double b, double b_total) {
                                                        double m_val = m_total <= 0 ? .5 : m / m_total;
                                                        double b_val = b_total <= 0 ? .5 : b / b_total;
                                                        return 1.0 - std::abs(.5 - m_val) - std::abs(.5 - b_val);
                                                    });
                halfspace2_t split_plane = plane;
                (void) hTotal;
                auto add_method = [&](wpoint_list_t const &cell) {
                    if (!cell.empty()) {
                        split_node_t cell_n;
                        for (auto &p : cell) {
                            if (split_plane.contains(p)) {
                                std::get<1>(cell_n).emplace_back(p);
                            } else {
                                std::get<0>(cell_n).emplace_back(p);
                            }
                        }
                        split_queue.push_back(cell_n);
                        std::push_heap(split_queue.begin(), split_queue.end(), cmp);
                    } else {
//                        std::cout << upper.size() << std::endl;
//                        std::cout << lower.size() << std::endl;
//                        std::cout << std::endl;
                    }
                };
                add_method(upper);
                add_method(lower);
            }
        }

        wpoint_list_t sample;
        auto sample_cell = [& sample, &gen](wpoint_list_t const& cell) {
            if (cell.empty()){
                return;
            }
            double u_w = computeTotal(cell);
            std::uniform_int_distribution<size_t> dist_u(0, cell.size() - 1);
            auto pt = cell[dist_u(gen)];
            sample.emplace_back(u_w, pt[0], pt[1], pt[2]);
        };

        for (auto& node : split_queue) {
            auto& [u, l] = node;
            sample_cell(u);
            sample_cell(l);
        }
        return sample;
    }

//    def approxHamSandwitch(pt_sets, eps):
//        """
//        Takes a set of sets of points in 2 dimensions and
//                attempts to find a polynomial that splits eacfh one approximately.
//        Approximately in this case is dictated by the parameter eps.
//        So an eps=1/8 set would find a plane that splits each set between 7/16, 8/16. An eps = 1/4
//        will find an approximate split 3/8 and 5/8.
//        """
//    #degree = int(.5 + math.ceil(1 / 2.0 * (math.sqrt(len(pt_sets) * 8 + 9) - 3)))
//        degree = len(pt_sets) - 1
//    #degree = int(degree + .5)
//        while True:
//        pts = [pt_set[random.randint(0, pt_set.shape[0] - 1), :] for pt_set in pt_sets]
//        coef = poly2dFit(np.vstack(pts), degree)
//        new_pt_sets = []
//        bad_pt_set = []
//        for pt_set in pt_sets:
//        polynomial_values = poly2dEval(pt_set, coef, degree)
//        mask = polynomial_values >= 0.0
//        eps_frac = float(mask.sum()) / len(pt_set)
//
//        if eps > abs(1/2.0 - eps_frac):
//        new_pt_sets.append(pt_set[mask, :])
//        new_pt_sets.append(pt_set[np.logical_not(mask), :])
//        else:
//        bad_pt_set.append(pt_set)
//        if len(new_pt_sets) > len(pt_sets): # disected more than half
//        return coef, new_pt_sets, bad_pt_set



//    wpoint_list_t polynomial_sample(const wpoint_list_t &pts, size_t s_size) {
//        /*
//         * Constructs a polynomial partitioning of the points pts.
//         */
//        while (1) {
//
//        }
//    }

    //lpoint_list_t ham_tree_sample(const lpoint_list_t &pts, size_t s_size);
}
