//
// Created by mmath on 5/28/17.
//
#include <functional>
#include <tuple>
#include <unordered_set>
//#include <zlib.h>
#include <memory>

#include "SparseGrid.hpp"
#include "Range.hpp"
#include "RectangleScan.hpp"
#include "FunctionApprox.hpp"

namespace pyscan {


    template<typename F>
    void quantiles(wpoint_list_t const& pts, std::vector<double>& indices, size_t r, F f) {
        auto order = util::sort_permutation(pts.begin(), pts.end(), [&](pt2_t const& p1, pt2_t const& p2){
            return f(p1) < f(p2);
        });
        double total_weight = std::accumulate(pts.begin(), pts.end(), 0, [](double w, wpt2_t const& p){
            return w + p.get_weight();
        });

        double eps_s = total_weight / r;
        double curr_weight = 0;
        std::for_each(order.begin(), order.end(), [&](size_t ix) {
            curr_weight += pts[ix].get_weight();
            if (curr_weight > eps_s) {
                indices.push_back(f(pts[ix]));
                curr_weight = 0;
            }
        });
    }

    Grid::Grid(wpoint_list_t const& red_points,
               wpoint_list_t const& blue_points) :
            r(red_points.size() + blue_points.size()),
            red_counts(r * r, 0),
            blue_counts(r * r, 0),
            x_coords(),
            y_coords() {

        //Compute the grid from the n_begin and n_end and then fill in the values with the two sampled things.
        auto getX = [](Point<> const &p1) {
            return p1(0);
        };
        auto getY = [](Point<> const &p1) {
            return p1(1);
        };

        for (auto & p : red_points){
            x_coords.push_back(getX(p));
            y_coords.push_back(getY(p));
        }
        for (auto & p : blue_points){
            x_coords.push_back(getX(p));
            y_coords.push_back(getY(p));
        }

        std::sort(x_coords.begin(), x_coords.end());
        std::sort(y_coords.begin(), y_coords.end());

        for (auto& point : red_points)  {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), point(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), point(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                red_counts[iy * r + ix] += point.get_weight();
            }
            total_red_weight += point.get_weight();
        }
        for (auto& point : blue_points) {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), point(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), point(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                blue_counts[iy * r + ix] += point.get_weight();
            }
            total_blue_weight += point.get_weight();
        }
    }

    /*
     * Takes a grid resolution,r, and a red set of points and a blue set of points. Computes a grid from this.z
     */
    Grid::Grid(size_t r_arg,
        wpoint_list_t const& red_points,
        wpoint_list_t const& blue_points) :
        r(r_arg),
        red_counts(r_arg * r_arg, 0),
        blue_counts(r_arg * r_arg, 0),
        x_coords(),
        y_coords()
    {

        //Compute the grid from the n_begin and n_end and then fill in the values with the two sampled things.
        auto getX = [](Point<> const &p1) {
            return p1(0);
        };
        auto getY = [](Point<> const &p1) {
            return p1(1);
        };



        quantiles(red_points, x_coords, r_arg, getX);
        quantiles(blue_points, x_coords, r_arg, getX);
        quantiles(red_points, y_coords, r_arg, getY);
        quantiles(blue_points, y_coords, r_arg, getY);
        std::sort(x_coords.begin(), x_coords.end());
        std::sort(y_coords.begin(), y_coords.end());

        for (auto& point : red_points)  {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), point(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), point(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                red_counts[iy * r + ix] += point.get_weight();
            }
            total_red_weight += point.get_weight();
        }
        for (auto& point : blue_points) {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), point(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), point(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                blue_counts[iy * r + ix] += point.get_weight();
            }
            total_blue_weight += point.get_weight();
        }
    }

    Grid::Grid(point_list_t const& net, wpoint_list_t const& red, wpoint_list_t const& blue) :
            r(net.size()),
            red_counts(r * r, 0),
            blue_counts(r * r, 0),
            x_coords(),
            y_coords() {

        for_each(net.begin(), net.end(), [&](pt2_t const& pt) {
            x_coords.push_back((pt)(0));
        });
        for_each(net.begin(), net.end(), [&] (pt2_t const& pt) {
            y_coords.push_back((pt)(1));
        });
        std::sort(x_coords.begin(), x_coords.end());
        std::sort(y_coords.begin(), y_coords.end());
        for (auto const& pt : red) {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), pt(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), pt(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                red_counts[iy * r + ix] += pt.get_weight();
            }
            total_red_weight += pt.get_weight();
        }
        for (auto const& pt : blue) {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), pt(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), pt(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                blue_counts[iy * r + ix] += pt.get_weight();
            }
            total_blue_weight += pt.get_weight();
        }
    }

    double Grid::totalRedWeight() const {
        return total_red_weight;
    }

    double Grid::totalBlueWeight() const {
        return total_blue_weight;
    }

    double Grid::redCount(size_t row, size_t col) const {
        return red_counts[row * r + col];
    }

    double Grid::blueCount(size_t row, size_t col) const {
        return blue_counts[row * r + col];
    }

    double Grid::redWeight(size_t row, size_t col) const {
        return red_counts[row * r + col] / (double) total_red_weight;
    }

    double Grid::blueWeight(size_t row, size_t col) const {
        return blue_counts[row * r + col] / (double) total_blue_weight;
    }

    double Grid::redSubWeight(Subgrid const &sg) const {
        double rc = 0;
        for (size_t i = sg.lowY(); i <= sg.upY(); i++) {
            for (size_t j = sg.lowX(); j <= sg.upX(); j++) {
                rc += redCount(i, j);
            }
        }
        return rc / (double) total_red_weight;
    }

    double Grid::blueSubWeight(Subgrid const &sg) const {
        double bc = 0;
        for (size_t i = sg.lowY(); i <= sg.upY(); i++) {
            for (size_t j = sg.lowX(); j <= sg.upX(); j++) {
                bc += blueCount(i, j);
            }
        }
        return bc / (double) total_blue_weight;
    }

    double Grid::yCoord(size_t row) const {
        return y_coords[row];
    }

    double Grid::xCoord(size_t col) const {
        return x_coords[col];
    }

    size_t Grid::size() const {
        return r;
    }

    Rectangle Grid::toRectangle(Subgrid const &sg) const {
        return Rectangle(xCoord(sg.upX()), yCoord(sg.upY()),
                         xCoord(sg.lowX()), yCoord(sg.lowY()));
    }

    /*
     * Simple 1/eps^4 algorithm described in the paper that just computes every subgrid.
     * This will work on a nonlinear function.
     */
    Subgrid max_subgrid(Grid const &grid, const discrepancy_func_t &func) {
        std::vector<double> red_count(grid.size(), 0);
        std::vector<double> blue_count(grid.size(), 0);


        auto t_red = static_cast<double>(grid.totalRedWeight());
        auto t_blue = static_cast<double>(grid.totalBlueWeight());

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
                    double red_l_count = 0;
                    double blue_l_count = 0;
                    for (size_t l = k; l < grid.size(); l++) {
                        red_l_count += red_count[l];
                        blue_l_count += blue_count[l];
                        double maxf = func(red_l_count,  t_red, blue_l_count, t_blue);
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


    Subgrid max_subgrid_linear(Grid const &grid, double a, double b) {
        std::vector<double> weight(grid.size(), 0);
        Subgrid max = Subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
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

    Subgrid max_subgrid_convex(Grid const &grid, double eps, discrepancy_func_t const &f) {
        /*
         * This uses the
         */
        auto phi = [&] (Vec2 const& v) {return f(v[0], 1.0,  v[1], 1.0); };
        Subgrid curr_subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
        Subgrid max_subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
        double maxV = 0;
        auto linemaxF = [&] (Vec2 const& dir) {
            curr_subgrid = max_subgrid_linear(grid, dir[0], dir[1]);
            Vec2 curr_mb{grid.redSubWeight(curr_subgrid), grid.blueSubWeight(curr_subgrid)};

            double curr_val = phi(curr_mb);
            if (curr_val > maxV) {
                max_subgrid = curr_subgrid;
                maxV = curr_val;
            }
            return curr_mb;
        };
        approximateHull(eps, phi, linemaxF);
        return max_subgrid;
    }

    size_t log2D (size_t val) {
        if (val == 1) return 0;
        unsigned int ret = 0;
        while (val > 1) {
            val >>= 1;
            ret++;
        }
        return ret;
    }

    size_t r_ratio(size_t n) {
        return n * log2D(n);
    }



    //////////////////////////////////////////////////////////////////////////////////
    //Begining OF LABELED RECTANGLE CODE//////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////




    class LabeledGrid {

        struct LabelW {
            size_t label;
            double weight;

            LabelW(size_t l, double w) : label(l), weight(w) {}
            LabelW(): label(0), weight(0) {}


            friend std::ostream& operator<<(std::ostream& os, const LabelW& el) {
                os << "(" << el.label << ", " << el.weight << ")";
                return os;
            }
        };


        using part_list_t = std::vector<double>;
        using grid_cell_t = std::vector<LabelW>;
        using grid_t = std::vector<grid_cell_t>;

        part_list_t x_values;
        part_list_t y_values;
        grid_t m_labels;
        grid_t b_labels;

    public:

        LabeledGrid(size_t r, lpoint_list_t m_points, lpoint_list_t b_points) {
            double m_total = computeTotal(m_points);
            double b_total = computeTotal(b_points);

            auto accum_and_partition = [&] (lpoint_list_t& labeled_points, double max_weight, size_t dim){

                /*
                 * I am going to assume that all the points have weight less than max_weight and then the
                 * weights are between max_weight and 2 * max_weight.
                 */
                std::sort(labeled_points.begin(), labeled_points.end(), [dim](Point<> const &p1, Point<> const &p2) {
                    return p1(dim) < p2(dim);
                });

                double curr_weight = 0;
                part_list_t partitions;
                std::unordered_set<size_t> active_label_set;
                for (auto& p : labeled_points) {
                    if (active_label_set.find(p.get_label()) == active_label_set.end()) {
                        curr_weight += p.get_weight();
                        active_label_set.emplace(p.get_label());
                    }
                    if (curr_weight > max_weight) {
                        partitions.emplace_back(p(dim));
                        curr_weight = 0;
                        active_label_set.clear();
                    }
                }
                return partitions;
            };

            auto update_parts = [&] (part_list_t &output, size_t dim) {
                auto part_m_dim = accum_and_partition(m_points, m_total / r, dim);
                auto part_b_dim = accum_and_partition(b_points, b_total / r, dim);
                output.resize(part_m_dim.size() + part_b_dim.size() + 2);
                *(output.begin()) = -std::numeric_limits<double>::infinity();
                std::merge(part_m_dim.begin(), part_m_dim.end(), part_b_dim.begin(), part_b_dim.end(), output.begin() + 1);
                *(output.end() - 1) = std::numeric_limits<double>::infinity();
            };

            update_parts(x_values, 0);
            update_parts(y_values, 1);

            m_labels.resize((x_values.size() - 1) * (y_values.size() - 1));
            b_labels.resize((x_values.size() - 1) * (y_values.size() - 1));

            auto grid_insert = [&] (grid_t &cells, lpoint_list_t const& lbl_pts) {
                for (auto &pt : lbl_pts) {
                    long ix = std::upper_bound(x_values.begin(), x_values.end(), pt(0)) - x_values.begin() - 1;
                    long iy = std::upper_bound(y_values.begin(), y_values.end(), pt(1)) - y_values.begin() - 1;
                    if (iy == (long)y_values.size() - 1 || ix == (long)x_values.size() - 1) {
                        //This can only happen if pt0 or pt1 is inf so we just clamp it back to fit in the box.
                        ix = ix - 1;
                        iy = iy - 1;
                    }
                    cells[iy * (x_values.size() - 1) + ix].emplace_back(pt.get_label(), pt.get_weight());
                }

            };
            grid_insert(m_labels, m_points);
            grid_insert(b_labels, b_points);
        }

        grid_cell_t const& get_m(size_t ix, size_t iy) const {
            return m_labels[iy * (x_values.size() - 1) + ix];
        }

        grid_cell_t const& get_b(size_t ix, size_t iy) const {
            return b_labels[iy * (x_values.size() - 1) + ix];
        }

        double x_val(size_t ix) const {
            return x_values[ix];
        }

        double y_val(size_t iy) const {
            return y_values[iy];
        }

        size_t x_size() const {
            return x_values.size() - 1;
        }

        size_t y_size() const {
            return y_values.size() - 1;
        }
    };

    std::tuple<Rectangle, double> max_rect_labeled(size_t r, double max_w,
                                                   lpoint_list_t const& m_points,
                                                   lpoint_list_t const& b_points,
                                                   double m_Total,
                                                   double b_Total,
                                                   const discrepancy_func_t& func) {

        LabeledGrid grid(r, m_points, b_points);

        Rectangle maxRect(0.0, 0.0, 0.0, 0.0);
        double max_stat = 0;

        for (size_t lower_j = 0; lower_j < grid.y_size(); lower_j++) {
            std::vector<std::unordered_map<size_t, double>> m_columns(grid.x_size(), std::unordered_map<size_t, double>());
            std::vector<std::unordered_map<size_t, double>> b_columns(grid.x_size(), std::unordered_map<size_t, double>());
            for (size_t upper_j = lower_j; upper_j < grid.y_size() &&
                                           std::abs(grid.y_val(upper_j + 1) - grid.y_val(lower_j)) < max_w; upper_j++) {
                for (size_t left_i = 0; left_i < grid.x_size(); left_i++) {

                    for (auto& lw : grid.get_m(left_i, upper_j)) {
                        m_columns[left_i][lw.label] = lw.weight;
                    }
                    for (auto& lw : grid.get_b(left_i, upper_j)) {
                        b_columns[left_i][lw.label] = lw.weight;
                    }
                }

                for (size_t left_i = 0; left_i < grid.x_size(); left_i++) {
                    //make sweep
                    std::unordered_set<size_t> curr_m_labels;
                    std::unordered_set<size_t> curr_b_labels;
                    double m_weight = 0;
                    double b_weight = 0;
                    for (size_t right_i = left_i; right_i < grid.x_size() &&
                                                  std::abs(grid.x_val(right_i + 1) - grid.x_val(left_i)) < max_w; right_i++) {
                        for (auto& label : m_columns[right_i]) {
                            if (curr_m_labels.find(label.first) == curr_m_labels.end()) {
                                m_weight += label.second;
                                curr_m_labels.emplace(label.first);
                            }
                        }
                        for (auto& label : b_columns[right_i]) {
                            if (curr_b_labels.find(label.first) == curr_b_labels.end()) {
                                b_weight += label.second;
                                curr_b_labels.emplace(label.first);
                            }
                        }
//                        Rectangle rect(grid.x_val(right_i + 1), grid.y_val(upper_j + 1), grid.x_val(left_i), grid.y_val(lower_j));
//                        double m_tmp = range_weight(rect, m_points);
//                        double b_tmp = range_weight(rect, b_points);
//                        std::cout << m_weight << " " << m_tmp << std::endl;
//                        std::cout << b_weight << " " << b_tmp << std::endl;
//                        std::cout << std::endl;

                        if (func(m_weight, m_Total, b_weight, b_Total) > max_stat) {
                            max_stat = func(m_weight, m_Total, b_weight, b_Total);
                            maxRect = Rectangle(grid.x_val(right_i + 1), grid.y_val(upper_j + 1), grid.x_val(left_i), grid.y_val(lower_j));
                        }
                    }
                }
            }
        }
        return std::make_tuple(maxRect, max_stat);
    }


    std::tuple<Rectangle, double> max_rect_labeled(size_t r, double max_w,
                                                   lpoint_list_t const& m_points,
                                                   lpoint_list_t const& b_points,
                                                   const discrepancy_func_t& func) {

        double m_Total = computeTotal(m_points);
        double b_Total = computeTotal(b_points);
        return max_rect_labeled(r, max_w, m_points, b_points, m_Total, b_Total, func);
    }


    std::tuple<Rectangle, double> max_rect_labeled_scale(
            size_t r,
            double max_r,
            double alpha,
            const point_list_t &net,
            const lpoint_list_t &red,
            const lpoint_list_t &blue,
            const discrepancy_func_t &f) {


        double red_tot = computeTotal(red);
        double blue_tot = computeTotal(blue);
        size_t grid_r = lround(floor(1 / alpha));
        SparseGrid<pt2_t> grid_net(net, grid_r);
        SparseGrid<lpt2_t> grid_red(red, grid_r), grid_blue(blue, grid_r);
        size_t sub_grid_size = lround(ceil(alpha / max_r));
        Rectangle max_rect;
        double max_stat = 0.0;
        for (auto center_cell = grid_net.begin(); center_cell != grid_net.end();) {
            std::vector<lpt2_t> net_chunk;
            std::vector<lpt2_t> red_chunk;
            std::vector<lpt2_t> blue_chunk;
            net_chunk.clear();
            red_chunk.clear();
            blue_chunk.clear();
            size_t i, j;
            std::tie(i, j) = grid_net.get_cell(center_cell->second);

            size_t end_k = i + sub_grid_size < grid_r ? i + sub_grid_size : grid_r;
            size_t end_l = j + sub_grid_size < grid_r ? j + sub_grid_size : grid_r;

            for (size_t k = i; k <= end_k; ++k) {
                for (size_t l = j; l <= end_l; ++l) {
                    auto red_range = grid_red(k, l);
                    for (auto it = red_range.first; it != red_range.second; ++it)
                        red_chunk.emplace_back(it->second);


                    auto blue_range = grid_blue(k, l);
                    for (auto it = blue_range.first; it != blue_range.second; ++it)
                        blue_chunk.emplace_back(it->second);
                }
            }
            auto [new_rect, local_max_stat] = max_rect_labeled(r, max_r, red_chunk, blue_chunk, red_tot, blue_tot, f);
            if (local_max_stat > max_stat) {
                max_rect = new_rect;
                max_stat = local_max_stat;
            }

            auto last = center_cell->first;
            do {
                ++center_cell;
            } while (center_cell != grid_net.end() && center_cell->first == last);
        }
        return std::make_tuple(max_rect, max_stat);
    }


    //////////////////////////////////////////////////////////////////////////////////
    //END OF LABELED RECTANGLE CODE//////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////



    //So each slab has list of splits and approximation list.
    //--splits are pointers to points.
    //--Approximation lists are lists of pointers as well..

    //Idea 2 use pointers. (so list of pointers to each point.) So at each slab have split list of pointers and
    // a merge list of pointers.

     class Interval {
         pt2_t* left;
         pt2_t* right;
         double value;
     public:
         Interval() : left(nullptr), right(nullptr), value(0.0) {}
         Interval(pt2_t* l, pt2_t* r, double v) : left(l), right(r), value(v) {}
         Interval(pt2_t* pt, double weight) : left(pt), right(pt), value(weight) {}

         pt2_t* get_r() const { return right; }
         pt2_t* get_l() const { return left; }
         double get_v() const {return value; }
     };

     class MaxIntervalAlt {
         Interval left_max;
         Interval right_max;
         Interval center_max;
     public:

         MaxIntervalAlt(pt2_t* pt, double weight) : left_max(pt, weight), right_max(pt, weight), center_max(pt, weight) {}

         MaxIntervalAlt(pt2_t* lpt, pt2_t* rpt) : left_max(lpt, lpt, 0.0), right_max(rpt, rpt, 0.0), center_max(lpt, lpt, 0.0) {}

         MaxIntervalAlt &operator+=(const MaxIntervalAlt& op) {
             if (center_max.get_v() < op.center_max.get_v()) {
                 center_max = op.center_max;
             }

             if (right_max.get_v() + op.left_max.get_v() > center_max.get_v()) {
                 center_max = Interval(right_max.get_l(), op.left_max.get_r(), right_max.get_v() + op.left_max.get_v());
             }

             //If left max is whole interval then we might need to extend the interval.
             if (left_max.get_r() == right_max.get_r()) {
                if (left_max.get_v() + op.left_max.get_v() > left_max.get_v()) {
                    left_max = Interval(left_max.get_l(), op.left_max.get_r(), left_max.get_v() + op.left_max.get_v());
                }
             }

             // op.right max is whole interval of op then we might need to extend op.right_max to include right_max
             if (op.right_max.get_l() == op.left_max.get_l() &&
                    op.right_max.get_v() + right_max.get_v() > op.right_max.get_v()) {
                 right_max = Interval(right_max.get_l(), op.right_max.get_r(), right_max.get_v() + op.right_max.get_v());
             } else {
                 right_max = op.right_max;
             }
             return *this;
         }

         pt2_t* left() const {
             return left_max.get_l();
         }
         pt2_t* right() const {
             return left_max.get_r();
         }

         Interval get_max() const {
             return center_max;
         }

     };


    merge_list_t compress(wpoint_it_t b, wpoint_it_t e, double m_w) {
        /*
         * Takes a vector of weighted points and constructed a smaller set of weighted points equally spaced with
         * the a different set of weights.
         */

//        if (wpts.empty()) {
//            return std::make_tuple(std::vector<pt2_t*>(), std::vector<double>());
//        }
//        double curr_weight = (*wpts.begin())->get_weight();
//        std::vector<double> weights;
//        std::vector<pt2_t *> break_pts;
//        for (auto it = wpts.begin() + 1; it != wpts.end(); ++it) {
//            if ((curr_weight + (*it)->get_weight()) > m_w) {
//                weights.emplace_back(curr_weight);
//                break_pts.emplace_back(*(it - 1));
//                curr_weight = 0;
//            }
//            curr_weight += (*it)->get_weight();
//        }
//        if (curr_weight == m_w) {
//            weights.emplace_back(curr_weight);
//            break_pts.emplace_back(*(wpts.end() - 1));
//        }
//        return std::make_tuple(break_pts, weights);
        merge_list_t merges;
        for (; b != e; ++b) {
            merges.emplace_back((*b)(0), 1.0);
        }
        return merges;
    }

    void normalize(wpoint_list_t& pts) {
        //Compute the total weight.
        double total_w = 0;
        for (auto& p : pts)  total_w += p.get_weight();
        std::transform(pts.begin(), pts.end(), pts.begin(), [total_w](const wpt2_t& pt){
            return wpt2_t(pt.get_weight() / total_w, pt[0], pt[1], pt[2]);
        });
    };


    Slab::Slab(std::weak_ptr<Slab> p,
            merge_list_t m_m,
            merge_list_t b_m,
            double ty,
            double by) :
                m_merges(std::move(m_m)),
                b_merges(std::move(b_m)),
                top_y(ty),
                bottom_y(by),
                parent(std::move(p)) {}


    double Slab::measure_interval(double mxx, double mnx, double a, double b) const {
        assert(mxx >= mnx);
        auto cmp = [](const merge_pair& p1, const merge_pair& p2) {
            return std::get<0>(p1) < std::get<0>(p2);
        };

        auto accum = [](double curr_val, const merge_pair& p2) {
            return curr_val + std::get<1>(p2);
        };
        auto m_u_b = std::lower_bound(m_merges.begin(), m_merges.end(), merge_pair(mxx, 0.0), cmp);
        auto m_l_b = std::lower_bound(m_merges.begin(), m_merges.end(), merge_pair(mnx, 0.0), cmp);
        double mval = a * std::accumulate(m_l_b, m_u_b, 0.0, accum);
        auto b_u_b = std::lower_bound(b_merges.begin(), b_merges.end(), merge_pair(mxx, 0.0), cmp);
        auto b_l_b = std::lower_bound(b_merges.begin(), b_merges.end(), merge_pair(mnx, 0.0), cmp);
        double bval = b * std::accumulate(b_l_b, b_u_b, 0.0, accum);
        return mval + bval;
    }

    SlabTree::SlabTree(std::vector<double> const &vert_decomp, wpoint_list_t ms, wpoint_list_t bs, double max_w) : mpts(
            std::move(ms)), bpts(std::move(bs)) {
        normalize(mpts);
        normalize(bpts);

        using wp_it = wpoint_list_t::iterator;
        using wrng_it = std::tuple<wp_it, wp_it>;
        using vert_it = std::vector<double>::const_iterator;
        using vrng_it = std::tuple<vert_it, vert_it>;

        using cell_t = std::tuple<slab_ptr, slab_ptr *, vrng_it, wrng_it, wrng_it>;
        using cell_list_t = std::vector<cell_t>;

        //Sort samples horizontally
        auto x_order = [](const pt2_t& p1, const pt2_t& p2) {
            return p1(0) < p2(0);
        };


        {
            double v_b = *(vert_decomp.begin());
            double v_e = *(vert_decomp.end() - 1);
            auto in_interval = [v_b, v_e] (const pt2_t& pt) {
                return v_b <= pt(1) && pt(1) < v_e;
            };
            auto m_new_end = std::partition(mpts.begin(), mpts.end(), in_interval);
            mpts.erase(m_new_end, mpts.end());
            auto b_new_end = std::partition(bpts.begin(), bpts.end(), in_interval);
            bpts.erase(b_new_end, bpts.end());
        }

        std::sort(mpts.begin(), mpts.end(), x_order);
        std::sort(bpts.begin(), bpts.end(), x_order);

        cell_list_t active;
        active.emplace_back(nullptr, &root,
                            std::make_tuple(vert_decomp.cbegin(), vert_decomp.cend()),
                            std::make_tuple(mpts.begin(), mpts.end()),
                            std::make_tuple(bpts.begin(), bpts.end()));

        std::vector<slab_ptr> roots;
        while (!active.empty()) {
            auto[parent, el_ptr, vrng, mrng, brng] = active.back();

            active.pop_back();

            auto[v_b, v_e] = vrng;
            auto[m_b, m_e] = mrng;
            auto[b_b, b_e] = brng;


            //Compress the current m_pts and b_pts and compute the weights.
            auto m_merges = compress(m_b, m_e, max_w);
            auto b_merges = compress(b_b, b_e, max_w);

            //Compute the slab now.
            *el_ptr = std::shared_ptr<Slab>(new Slab(parent, m_merges, b_merges, *(v_e - 1), *v_b));

            //Now emplace back to the active
            if ((v_e - v_b) > 2) {

                auto v_mid = (v_b + (v_e - v_b) / 2);
                assert(v_mid != v_e);
                auto v_break = *v_mid;
                auto part_f = [v_break](const pt2_t& ptr) {
                    return ptr(1) < v_break;
                };
                auto m_iter_splt = std::stable_partition(m_b, m_e, part_f);
                auto b_iter_splt = std::stable_partition(b_b, b_e, part_f);
                active.emplace_back(*el_ptr,
                                    &((*el_ptr)->down),
                                    std::make_tuple(v_b, v_mid + 1),
                                    std::make_tuple(m_b, m_iter_splt),
                                    std::make_tuple(b_b, b_iter_splt));
                active.emplace_back(*el_ptr,
                                    &((*el_ptr)->up),
                                    std::make_tuple(v_mid, v_e),
                                    std::make_tuple(m_iter_splt, m_e),
                                    std::make_tuple(b_iter_splt, b_e));
            } else {
                roots.emplace_back(*el_ptr);
            }
        }

        auto order_f = [] (const merge_pair& mp1, const merge_pair& mp2) {
            return std::get<0>(mp1) < std::get<0>(mp2);
        };

        //Now go back through and create a list of all the lower splits.
        while (!roots.empty()) {
            auto child = roots.back();
            roots.pop_back();

            //First merge the m_merges and b_merges
            child->split_offsets.reserve(child->split_offsets.size() + child->m_merges.size() + child->b_merges.size());

            std::merge(child->m_merges.begin(),
                       child->m_merges.end(),
                       child->b_merges.begin(),
                       child->b_merges.end(),
                       std::back_inserter(child->split_offsets), order_f);
            auto mid = child->split_offsets.begin() + child->split_offsets.size();
            //Now merge into the already existing split offsets.
            std::inplace_merge(child->split_offsets.begin(), mid, child->split_offsets.end(), order_f);

            //Now merge into the parent.
            auto p = child->parent.lock();
            if (p != nullptr) {
                size_t offset = p->split_offsets.size();
                std::copy(child->split_offsets.begin(), child->split_offsets.end(),
                          std::back_inserter(p->split_offsets));
                std::inplace_merge(p->split_offsets.begin(), p->split_offsets.begin() + offset, p->split_offsets.end(), order_f);
            }
        }
    }

    SlabTree::slab_ptr SlabTree::get_containing(Rectangle const &rect) const {
        auto curr_root = root;
        while (curr_root != nullptr) {
            //Check if this region is contained completely in the left or right branch.
            if (curr_root->has_mid()) {
                double mid = curr_root->get_mid();
                if (rect.upY() < mid) {
                    curr_root = curr_root->down;
                } else if (rect.lowY() > mid) {
                    curr_root = curr_root->up;
                } else {
                    return curr_root;
                }
            } else {
                return curr_root;
            }
        }
        return curr_root;
    }

    double SlabTree::measure_rect(Rectangle const &rect, double a, double b) const {
        auto curr_root = get_containing(rect);

        if (curr_root == nullptr) {
            return 0.0;
        }
        auto& upper_bound = curr_root->up;
        double sum = 0;
        while (upper_bound != nullptr) {

            double midpoint;
            if (upper_bound->down != nullptr) midpoint = upper_bound->down->top_y;
            else if (upper_bound->up != nullptr) midpoint = upper_bound->up->bottom_y;
            else {
                //sum += upper_bound->measure_interval(rect.upX(), rect.lowX(), a, b);
                break;
            }

            if (midpoint < rect.upY()) {
                sum += !upper_bound->down ? 0.0 : upper_bound->down->measure_interval(rect.upX(), rect.lowX(), a, b);
                upper_bound = upper_bound->up;
            } else {
                upper_bound = upper_bound->down;
            }
        }
        auto& lower_bound = curr_root->down;
        while (lower_bound != nullptr) {
            double midpoint;
            if (lower_bound->down != nullptr) midpoint = lower_bound->down->top_y;
            else if (lower_bound->up != nullptr) midpoint = lower_bound->up->bottom_y;
            else {
                //sum += lower_bound->measure_interval(rect.upX(), rect.lowX(), a, b);
                break;
            }

            if (midpoint >= rect.lowY()) {
                sum += !lower_bound->up ? 0.0 : lower_bound->up->measure_interval(rect.upX(), rect.lowX(), a, b);
                lower_bound = lower_bound->down;
            } else {
                lower_bound = lower_bound->up;
            }
        }
        return sum;
    }


//    std::tuple<Rectangle, double> SlabTree::max_rectangle(double m_a, double b_b) {
//        //Initialize with list of maximum intervals.
//        assert(!slab_heap.empty());
//        auto &curr_splits = slab_heap[0].split_offsets;
//        assert(!curr_splits.empty());
//
//        std::vector<MaxIntervalAlt> initial_intervals;
//        initial_intervals.reserve(curr_splits.size() - 1);
//
//        for (size_t i = 0; i < curr_splits.size() - 1; ++i) {
//            initial_intervals.emplace_back(curr_splits[i], curr_splits[i + 1]);
//        }
//
//        std::vector<std::tuple<Slab *, Slab *, Slab *, Slab *, pt2_t *, pt2_t *, std::vector<MaxIntervalAlt>>> slab_stack;
//
//        assert(slab_heap[0].up != nullptr && slab_heap[0].down != nullptr);
//
//        slab_stack.emplace_back(slab_heap[0].up->up, slab_heap[0].up->down,
//                                slab_heap[0].down->up, slab_heap[0].down->down,
//                                slab_heap[0].up->top_y, slab_heap[0].down->bottom_y,
//                                initial_intervals
//        );
//
//
//        auto update_weight = [](std::vector<MaxIntervalAlt> &max_intervals, std::vector<pt2_t *> merges,
//                                std::vector<double> weights, double scale) {
//            size_t j = 0;
//            for (size_t i = 0; i < merges.size(); ++i) {
//                while (j < max_intervals.size() && merges[i] != max_intervals[j].right()) {
//                    j++;
//                }
//                max_intervals[j] += MaxIntervalAlt(merges[i], scale * weights[i]);
//            }
//        };
//
//        auto merge_splits = [](std::vector<MaxIntervalAlt> const &max_intervals, std::vector<pt2_t *> const &splits) {
//            size_t j = 0;
//            std::vector<MaxIntervalAlt> new_intervals;
//            for (size_t i = 0; i < splits.size(); ++i) {
//                while (j < max_intervals.size() - 1 && splits[i] != max_intervals[j].right()) {
//                    new_intervals.emplace_back(max_intervals[j]);
//                    j++;
//                }
//                new_intervals.emplace_back(max_intervals[j]);
//                new_intervals[j] += max_intervals[j + 1];
//                j++;
//            }
//            return new_intervals;
//        };
//
//        auto update_and_merge = [&](Slab *top, Slab *bottom, std::vector<MaxIntervalAlt> &m1) {
//            auto m2 = m1;
//            if (bottom != nullptr) { //if bottom is nullptr then we merge top or don't merge top.
//                update_weight(m1, bottom->m_merges, bottom->m_weights, m_a);
//                update_weight(m1, bottom->b_merges, bottom->b_weights, b_b);
//                m1 = merge_splits(m1, bottom->split_offsets);
//
//                m2 = merge_splits(m2, bottom->m_merges);
//                m2 = merge_splits(m2, bottom->b_merges);
//            }
//            if (top != nullptr) {
//                m1 = merge_splits(m1, top->m_merges);
//                m1 = merge_splits(m1, top->b_merges);
//                m2 = merge_splits(m2, top->split_offsets);
//            }
//            return m2;
//        };
//        auto get_top = [](Slab *p) {
//            return p == nullptr ? nullptr : p->up;
//        };
//        auto get_bottom = [](Slab *p) {
//            return p == nullptr ? nullptr : p->down;
//        };
//
//        auto get_mid = [](Slab *top, Slab *bottom) {
//            //This assumes that both bottom and top cannot be null.
//            return bottom == nullptr ? top->bottom_y : bottom->top_y;
//        };
//
//        double max_v = 0;
//        Rectangle max_rect;
//        while (!slab_stack.empty()) {
//            auto[top_top, top_bottom, bottom_top, bottom_bottom, p_up, p_low, max_intervals] = slab_stack.back();
//            slab_stack.pop_back();
//            if (max_intervals.size() > 1) {
//                // m4 doesn't have any merges
//                // m2 is bottom_top merged.
//                // m3 is top_bottom merged.
//                // max_intervals is top_bottom and bottom_top merged.
//                assert(!(top_top == nullptr &&
//                         top_bottom == nullptr &&
//                         bottom_top == nullptr &&
//                         bottom_bottom == nullptr));
//
//                if (!(top_top == nullptr && top_bottom == nullptr)) {
//                    auto m2 = update_and_merge(top_top, top_bottom, max_intervals);
//                    if (!(bottom_bottom == nullptr && bottom_top == nullptr)) {
//                        auto m4 = update_and_merge(bottom_bottom, bottom_top, m2);
//                        slab_stack.emplace_back(get_top(top_bottom), get_bottom(top_bottom),
//                                                get_top(bottom_top), get_bottom(bottom_top),
//                                                get_mid(top_top, top_bottom), get_mid(bottom_top, bottom_bottom), m4);
//                    }
//                    slab_stack.emplace_back(get_top(top_bottom), get_bottom(top_bottom),
//                                            get_top(bottom_bottom), get_bottom(bottom_bottom),
//                                            get_mid(top_top, top_bottom), p_low, m2);
//                }
//                if (!(bottom_bottom == nullptr && bottom_top == nullptr)) {
//                    auto m3 = update_and_merge(bottom_bottom, bottom_top, max_intervals);
//                    slab_stack.emplace_back(get_top(top_top), get_bottom(top_top),
//                                            get_top(bottom_top), get_bottom(bottom_top),
//                                            p_up, get_mid(bottom_top, bottom_bottom), m3);
//                }
//                slab_stack.emplace_back(get_top(top_top), get_bottom(top_top),
//                                        get_top(bottom_bottom), get_bottom(bottom_bottom),
//                                        p_up, p_low, max_intervals);
//
//            } else {
//                if (max_v < max_intervals[0].get_max().get_v()) {
//                    double lx = max_intervals[0].get_max().get_l()->operator()(0);
//                    double rx = max_intervals[0].get_max().get_r()->operator()(0);
//                    double uy = (*p_up)(1);
//                    double ly = (*p_low)(1);
//                    max_rect = Rectangle(rx, uy, lx, ly);
//                    max_v = max_intervals[0].get_max().get_v();
//                }
//            }
//        }
//        return std::make_tuple(max_rect, max_v);
//
//    }



//    std::tuple<Rectangle, double> max_rectangle(const wpoint_list_t& pts, double eps) {
//
//        return tree.max_rectangle(a, b);
//    }
}
