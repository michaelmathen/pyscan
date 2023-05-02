//
// Created by mmath on 5/28/17.
//
#include <functional>
#include <tuple>
#include <unordered_set>
//#include <zlib.h>
#include <memory>
#include <unordered_map>

#include "Gridding.hpp"
#include "Range.hpp"
#include "RectangleScan.hpp"
#include "FunctionApprox.hpp"

namespace pyscan {



    point_list_t quantiles(wpoint_list_t const& pts, double m_w) {
        double curr_w = 0;
        point_list_t break_pts;
        for (auto b = pts.begin(); b != pts.end(); ++b) {
            if ((curr_w + b->get_weight()) > m_w) {
                break_pts.emplace_back((*b)[0], (*b)[1], (*b)[2]);
                curr_w = 0;
            } else {
                curr_w += b->get_weight();
            }
        }
        return break_pts;
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
        wpoint_list_t red_points,
        wpoint_list_t blue_points) :
        r(r_arg),
        red_counts(r_arg * r_arg, 0),
        blue_counts(r_arg * r_arg, 0),
        x_coords(),
        y_coords(),
        total_red_weight(computeTotal(red_points)),
        total_blue_weight(computeTotal(blue_points))
    {

        //Compute the grid from the n_begin and n_end and then fill in the values with the two sampled things.
        auto xcmp = [](Point<> const &p1, Point<> const& p2) {
            return p1(0) < p2(0);
        };
        auto ycmp = [](Point<> const &p1, Point<> const& p2) {
            return p1(1) < p2(1);
        };


        auto insert_coords = [](point_list_t const& pts, size_t dim, std::vector<double>& coords) {
            coords.reserve(pts.size() + coords.size());
            for (auto& pt : pts) {
                coords.emplace_back(pt(dim));
            }
        };

        std::sort(red_points.begin(), red_points.end(), xcmp);
        std::sort(blue_points.begin(), blue_points.end(), xcmp);
        auto r_p_x = quantiles(red_points, 2 * total_red_weight / r_arg);
        auto b_p_x = quantiles(blue_points, 2 * total_blue_weight / r_arg);

        insert_coords(r_p_x, 0, x_coords);
        insert_coords(b_p_x, 0, x_coords);

        std::sort(red_points.begin(), red_points.end(), ycmp);
        std::sort(blue_points.begin(), blue_points.end(), ycmp);
        auto r_p_y = quantiles(red_points, 2 * total_red_weight / r_arg);
        auto b_p_y = quantiles(blue_points, 2 * total_blue_weight / r_arg);

        insert_coords(r_p_y, 1, y_coords);
        insert_coords(b_p_y, 1, y_coords);

        std::inplace_merge(x_coords.begin(), x_coords.begin() + r_p_x.size(), x_coords.end());
        std::inplace_merge(y_coords.begin(), y_coords.begin() + r_p_y.size(), y_coords.end());


        for (auto& point : red_points)  {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), point(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), point(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                red_counts[iy * r + ix] += point.get_weight();
            }
        }
        for (auto& point : blue_points) {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), point(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), point(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                blue_counts[iy * r + ix] += point.get_weight();
            }
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

        Rectangle max_rect;
        double max_stat = 0.0;
        auto bbop = bbox(net, red, blue);
        if (!bbop.has_value()) {
            return std::make_tuple(max_rect, max_stat);
        }
        auto bb = bbop.value();

        double red_tot = computeTotal(red);
        double blue_tot = computeTotal(blue);
        SparseGrid<pt2_t> grid_net(bb, net, alpha);
        SparseGrid<lpt2_t> grid_red(bb, red, alpha), grid_blue(bb, blue, alpha);
        auto grid_r = grid_net.get_grid_size();
        size_t sub_grid_size = lround(ceil(alpha / max_r));

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
    //Max Rectangle Code//////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////



    template<typename Arr>
    decltype(auto) access(Arr& arr1, Arr& arr2, size_t i) {
        if (arr1.size() <= i) {
            return arr2.at(i - arr1.size());
        } else {
            return arr1.at(i);
        }
    }

    std::tuple<epoint_list_t, epoint_list_t, std::unordered_map<size_t, double>, std::unordered_map<size_t, double> >
            to_epoints(wpoint_list_t const& mpts, wpoint_list_t const& bpts){
        /*
         * Maps every point to a new point with integer index and also creates a reverse mapping so
         * we can go from integer index back to the floating point value.
         */

        auto get_ix_list = [&] (size_t dim) {
            auto cmp = [&](size_t i1, size_t i2) {
                return access(mpts, bpts, i1)(dim) < access(mpts, bpts, i2)(dim);
            };
            std::vector<size_t> p(mpts.size() + bpts.size());
            std::iota(p.begin(), p.end(), 0);
            std::sort(p.begin(), p.end(), cmp);
            return p;
        };
        auto ixx = get_ix_list(0);
        auto ixy = get_ix_list(1);
        epoint_list_t empts(mpts.size());
        epoint_list_t ebpts(bpts.size());

        std::unordered_map<size_t, double> x_map;
        std::unordered_map<size_t, double> y_map;
        for (size_t i = 0; i < mpts.size() + bpts.size(); i++) {
            auto& tmp = access(empts, ebpts, ixx[i]);
            tmp(0) = i + 1;

            x_map[i + 1] = access(mpts, bpts, ixx[i])(0);
            y_map[i + 1] = access(mpts, bpts, ixy[i])(1);

            tmp.set_weight(access(mpts, bpts, ixx[i]).get_weight());
            access(empts, ebpts, ixy[i])(1) = i + 1;
        }
        return std::make_tuple(empts, ebpts, x_map, y_map);
    }




     epoint_list_t even_compress_internal(epoint_it_t b, epoint_it_t e, double m_w) {
        /*
         * Takes a vector of weighted points and constructed a smaller set of weighted points equally spaced with
         * the a different set of weights.
         */
        double curr_w = 0;
        std::vector<ept_t> break_pts;
        for (; b != e; ++b) {
            if ((curr_w + b->get_weight()) > m_w) {
                break_pts.emplace_back(b->get_x(), b->get_y(), curr_w + b->get_weight());
                curr_w = 0;
            } else {
                curr_w += b->get_weight();
            }
        }
        return break_pts;
    }


    template <typename UURG>
    epoint_list_t block_compress_internal(epoint_it_t b, epoint_it_t e, double m_w, UURG&& gen) {
        /*
         * Takes a vector of weighted points and constructed a smaller set of weighted points equally spaced with
         * the a different set of weights.
         */
        double curr_w = 0;
        std::vector<ept_t> break_pts;
        auto last_b = b;
        for (; b != e; ++b) {
            if ((curr_w + b->get_weight()) > m_w) {
                std::uniform_int_distribution<size_t> dist(0, static_cast<size_t>(b - last_b ));
                break_pts.emplace_back(*(last_b + dist(gen)));
                break_pts.back().set_weight(curr_w + b->get_weight());
                curr_w = 0;
            } else {
                curr_w += b->get_weight();
            }
        }
        return break_pts;
    }

//    epoint_list_t cascade_compress_internal(
//            epoint_it_t b,
//            epoint_it_t e,
//            epoint_it_t uncle_comp_b,
//            epoint_it_t uncle_comp_e,
//            epoint_it_t uncle_b,
//            epoint_it_t uncle_e,
//            double m_w) {
//        for (auto b1 = uncle_comp_b; b1 != uncle_comp_e; ++b1) {
//            b1->set_weight(-b1->get_weight());
//        }
//        epoint_list_t compressed;
//        compressed.resize((e - b) + (uncle_comp_e - uncle_comp_b) + (uncle_e - uncle_b));
//        std::merge(b, e, uncle_comp_b, uncle_comp_e, std::back_inserter(compressed));
//        std::copy(uncle_b, uncle_e, std::back_inserter(compressed));
//        std::inplace_merge(compressed.begin(), compressed.begin() + (e - b) + (uncle_comp_e - uncle_comp_b), compressed.end());
//        return even_compress_internal(compressed.begin(), compressed.end(), m_w);
//    }

    void normalize(wpoint_list_t& pts) {
        //Compute the total weight.
        double total_w = 0;
        for (auto& p : pts) {
            total_w += p.get_weight();
        }
        for (size_t i = 0; i < pts.size(); i++) {
            auto pt = pts[i];
            pts[i] = wpt2_t(pt.get_weight() / total_w, pt[0], pt[1], pt[2]);
        }
    };


    Slab::Slab(std::weak_ptr<Slab> p,
            epoint_list_t m_m,
            epoint_list_t b_m,
            size_t ty,
            size_t by) :
                m_merges(std::move(m_m)),
                b_merges(std::move(b_m)),
                top_y(ty),
                bottom_y(by),
                parent(std::move(p)) {}


    double Slab::measure_interval(size_t mxx, size_t mnx, double a, double b) const {
        /*
         * Measures the interval [mnx, mxx) associated with a single slab.
         */
        auto accum = [](double curr_val, const ept_t& p2) {
            return curr_val + p2.get_weight();
        };
        auto m_u_b = std::lower_bound(m_merges.begin(), m_merges.end(), EPoint(mxx, 0, 0.0));
        auto m_l_b = std::lower_bound(m_merges.begin(), m_merges.end(), EPoint(mnx, 0, 0.0));
        double mval = a * std::accumulate(m_l_b, m_u_b, 0.0, accum);
        auto b_u_b = std::lower_bound(b_merges.begin(), b_merges.end(), EPoint(mxx, 0, 0.0));
        auto b_l_b = std::lower_bound(b_merges.begin(), b_merges.end(), EPoint(mnx, 0, 0.0));
        double bval = b * std::accumulate(b_l_b, b_u_b, 0.0, accum);
        return mval + bval;
    }

    SlabTree::SlabTree(epoint_list_t mpts, epoint_list_t bpts, double max_w) :
        total_m(computeTotal(mpts)), total_b(computeTotal(bpts)) {
        /*
         * Create a vertical decomposition of the point set so that we can split the points into a sequence of horizontal
         * strips where each contains at most max_w.
         */
        auto y_order = [](const ept_t& p1, const ept_t& p2){
            return p1.get_y() < p2.get_y();
        };
        std::sort(mpts.begin(), mpts.end(), y_order);
        std::sort(bpts.begin(), bpts.end(), y_order);
        std::vector<ept_t> merge_buffer;
        merge_buffer.reserve(mpts.size() + bpts.size());
        std::merge(mpts.begin(), mpts.end(), bpts.begin(), bpts.end(), std::back_inserter(merge_buffer), y_order);

        std::vector<size_t> vert_decomp;
        vert_decomp.reserve((size_t) (1 / max_w));
        vert_decomp.emplace_back(0);
        double curr_weight = 0;
        for (auto &p : merge_buffer) {
            curr_weight += p.get_weight();
            if (curr_weight >= max_w * (total_m + total_b)) {
                curr_weight = 0;
                vert_decomp.emplace_back(p.get_y());
            }
        }
        auto mmx_it = mpts.end() - 1;
        auto bmx_it = bpts.end() - 1;
        if (mmx_it->get_y() < bmx_it->get_y()) {
            vert_decomp.emplace_back(bmx_it->get_y() + 1);
        } else {
            vert_decomp.emplace_back(mmx_it->get_y() + 1);
        }

        SlabTree::init(mpts, bpts, vert_decomp);
    }

    SlabTree::SlabTree(std::vector<size_t> vert_decomp, epoint_list_t mpts, epoint_list_t bpts) :
        total_m(computeTotal(mpts)), total_b(computeTotal(bpts)) {
        /*
         * Create a vertical decomposition of the point set so that we can split the points into a sequence of horizontal
         * strips where each contains at most max_w.
         */

        vert_decomp.insert(vert_decomp.begin(), 0);
        auto y_cmp = [] (ept_t const& p1, ept_t const& p2) {
            return p1.get_y() < p2.get_y();
        };
        auto mmx_it = max_element(mpts.begin(), mpts.end(), y_cmp);
        auto bmx_it = max_element(bpts.begin(), bpts.end(), y_cmp);
        if (mmx_it->get_y() < bmx_it->get_y()) {
            vert_decomp.emplace_back(bmx_it->get_y() + 1);
        } else {
            vert_decomp.emplace_back(mmx_it->get_y() + 1);
        }
        SlabTree::init(std::move(mpts), std::move(bpts), vert_decomp);
    }

    void SlabTree::init(epoint_list_t mpts, epoint_list_t bpts, std::vector<size_t> const& vert_decomp) {

        //std::sort(vert_decomp.begin(), vert_decomp.end());


        using wrng_it = std::tuple<epoint_it_t , epoint_it_t >;
        using vert_it = std::vector<size_t>::const_iterator;
        using vrng_it = std::tuple<vert_it, vert_it>;

        using cell_t = std::tuple<slab_ptr, slab_ptr *, vrng_it, wrng_it, wrng_it>;
        using cell_list_t = std::vector<cell_t>;


        //Sort these by x-axis
        std::sort(mpts.begin(), mpts.end());
        std::sort(bpts.begin(), bpts.end());

        cell_list_t active;
        active.emplace_back(nullptr, &root,
                            std::make_tuple(vert_decomp.cbegin(), vert_decomp.cend()),
                            std::make_tuple(mpts.begin(), mpts.end()),
                            std::make_tuple(bpts.begin(), bpts.end()));

        std::queue<slab_ptr> roots;
        while (!active.empty()) {
            auto[parent, el_ptr, vrng, mrng, brng] = active.back();

            active.pop_back();

            auto[v_b, v_e] = vrng;
            auto[m_b, m_e] = mrng;
            auto[b_b, b_e] = brng;


            //Compress the current m_pts and b_pts and compute the weights.
            //Need the
            auto m_merges = epoint_list_t(m_b, m_e);
            auto b_merges = epoint_list_t(b_b, b_e);

            //Compute the slab now.
            *el_ptr = std::shared_ptr<Slab>(new Slab(parent, m_merges, b_merges, *(v_e - 1), *v_b));

            //Now emplace back to the active
            if ((v_e - v_b) > 2) {

                auto v_mid = (v_b + (v_e - v_b) / 2);
                assert(v_mid != v_e);
                auto v_break = *v_mid;
                auto part_f = [v_break](const ept_t& ptr) {
                    //Has to be < so that it matches ERectangle contains
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
            }
        }
    }


    std::vector<slab_ptr> SlabTree::get_leaves() const {
        std::vector<decltype(root)> curr_stack;
        curr_stack.push_back(root);
        std::vector<decltype(root)> leaves;
        while (!curr_stack.empty()) {
            auto el = curr_stack.back();
            curr_stack.pop_back();
            if (el->up != nullptr) {
                curr_stack.push_back(el->up);
            }
            if (el->down != nullptr) {
                curr_stack.push_back(el->down);
            }
            if (!el->up  && !el->down){
                leaves.push_back(el);
            }
        }
        return leaves;
    }

    void SlabTree::reset_splits() {
        std::vector<decltype(root)> curr_stack;
        curr_stack.push_back(root);
        while (!curr_stack.empty()) {
            auto el = curr_stack.back();
            curr_stack.pop_back();
            if (el != nullptr) {
                el->global_split_offset = std::vector<size_t>();
                curr_stack.push_back(el->down);
                curr_stack.push_back(el->up);
            }
        }
    }

    void SlabTree::compute_splits() {

        auto roots = get_leaves();

        using s_it = std::vector<size_t>::iterator;
        auto inplace_set_union = [](std::vector<size_t>& memory, s_it b, s_it mid, s_it e) {
            // Now merge into the already existing split offsets.
            std::inplace_merge(b, mid, e);
            auto end_it = std::unique(b, e);
            memory.erase(end_it, e);
        };
        //Now go back through and create a list of all the lower splits.
        while (!roots.empty()) {
            auto child = roots.back();
            roots.pop_back();
            //Now merge into the parent.
            auto p = child->parent.lock();
            if (p != nullptr) {
                //First merge the m_merges and b_merges into the end of the existing global_split_offset
                std::vector<size_t> curr_splits;
                curr_splits.reserve(child->m_merges.size() + child->b_merges.size());
                for (auto& pt : child->m_merges) { curr_splits.emplace_back(pt.get_x()); }
                for (auto& pt : child->b_merges) { curr_splits.emplace_back(pt.get_x()); }

                inplace_set_union(curr_splits, curr_splits.begin(),
                                  curr_splits.begin() + child->m_merges.size(),
                                  curr_splits.end());

                size_t offset = curr_splits.size();
                curr_splits.reserve(curr_splits.size() + child->global_split_offset.size());
                std::copy(child->global_split_offset.begin(), child->global_split_offset.end(),
                          std::back_inserter(curr_splits));
                inplace_set_union(curr_splits, curr_splits.begin(), curr_splits.begin() + offset, curr_splits.end());

                size_t offset2 = p->global_split_offset.size();
                std::copy(curr_splits.begin(), curr_splits.end(), std::back_inserter(p->global_split_offset));

                inplace_set_union(p->global_split_offset, p->global_split_offset.begin(),
                                  p->global_split_offset.begin() + offset2, p->global_split_offset.end());
                if (p->up == nullptr || child == p->up) {
                    roots.emplace_back(p);
                }
            }
        }
    }

    void SlabTree::block_compress(double max_w) {

        std::vector<slab_ptr> curr_slabs;
        if (root == nullptr) {
            return;
        }
        std::random_device rd;
        std::minstd_rand gen(rd());

        curr_slabs.emplace_back(root->up);
        curr_slabs.emplace_back(root->down);
        while (!curr_slabs.empty()) {
            auto curr_root = curr_slabs.back();
            curr_slabs.pop_back();
            if (curr_root != nullptr) {
                curr_root->m_merges = block_compress_internal(curr_root->m_merges.begin(), curr_root->m_merges.end(), max_w * total_m, gen);
                curr_root->b_merges = block_compress_internal(curr_root->b_merges.begin(), curr_root->b_merges.end(), max_w * total_b, gen);
                curr_slabs.emplace_back(curr_root->up);
                curr_slabs.emplace_back(curr_root->down);
            }
        }
    }

    void SlabTree::even_compress(double max_w) {
        std::vector<slab_ptr> curr_slabs;
        curr_slabs.emplace_back(root);
        while (!curr_slabs.empty()) {
            auto curr_root = curr_slabs.back();
            curr_slabs.pop_back();
            if (curr_root != nullptr) {
                curr_root->m_merges = even_compress_internal(curr_root->m_merges.begin(), curr_root->m_merges.end(),
                                               max_w * total_m);
                curr_root->b_merges = even_compress_internal(curr_root->b_merges.begin(), curr_root->b_merges.end(), max_w * total_b);
                curr_slabs.emplace_back(curr_root->up);
                curr_slabs.emplace_back(curr_root->down);
            }
        }
    }

//    void SlabTree::cascade_compress(double eps){
//        /*
//         * This method will first partition m_merges and b_merges into log^4 ( 1/ eps ) partitions and sample
//         * a certain number of points from each cell. This will follow the scheme from A1 in
//         * Computing Approximate Statistical Discrepancy.
//         *
//         * Then we construct a 1/eps by 1/eps grid over the epoint m_merges and epoint b_merges.
//         * We then compute a compression for each slab while keeping track of the compressed representation
//         * of the points between the slab and the midline and a set of points corresponding to grid columns.
//         * We then use the scheme from A3 to create the compressed representation in each slab.
//         */
////        auto lg2 = [] (size_t val) {
////            int lg = 0;
////            while (val >>= 1) ++lg;
////            return lg;
////        };
////        auto special_lglg = [&lg2] (size_t size) {
////            size_t lg = std::max(lg2(size), 1);
////            //Get the correct power of c
////            return std::max(size / (lg * lg * lg * lg);
////        };
//
//
//
//    }

    slab_ptr SlabTree::get_containing(size_t upY, size_t lowY) const {
        assert(upY >= lowY);
        auto curr_root = root;
        while (curr_root != nullptr) {
            //Check if this region is contained completely in the left or right branch.
            if (curr_root->has_mid()) {
                size_t mid = curr_root->get_mid();

                if (upY < mid) {
                    curr_root = curr_root->down;
                } else if (lowY > mid) {
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

    std::vector<size_t> SlabTree::get_vert_breakpoints() const {
        std::vector<decltype(root)> curr_stack;
        curr_stack.push_back(root);
        std::vector<decltype(root)> leaves;
        while (!curr_stack.empty()) {
            auto el = curr_stack.back();
            curr_stack.pop_back();
            if (el->up != nullptr) {
                curr_stack.push_back(el->up);
            }
            if (el->down != nullptr) {
                curr_stack.push_back(el->down);
            }
            if (!el->up  && !el->down){
                leaves.push_back(el);
            }
        }

        std::vector<size_t> divisions;
        divisions.reserve(leaves.size() - 1);
        for (size_t i = 1; i < leaves.size(); i++) {
            divisions.push_back(leaves[i - 1]->top_y);
        }
        return divisions;
    }

    double SlabTree::measure_rect(ERectangle const &rect, double a, double b) const {
        auto curr_root = get_containing(rect.upY(), rect.lowY());
        if (curr_root == nullptr) {
            return 0.0;
        }
        auto upper_bound = curr_root->up;
        double rect_total = 0;
        while (upper_bound != nullptr) {
            if (upper_bound->has_mid()) {
                double midpoint = upper_bound->get_mid();
                if (midpoint < rect.upY()) {
                    double tmp = upper_bound->down == nullptr ? 0.0 : upper_bound->down->measure_interval(rect.upX(), rect.lowX(), a, b);
                    rect_total += tmp;
                    upper_bound = upper_bound->up;
                } else {
                    upper_bound = upper_bound->down;
                }
            } else {
                rect_total += upper_bound->top_y <= rect.upY() ? upper_bound->measure_interval(rect.upX(), rect.lowX(), a, b) : 0.0;
                break;
            }
        }
        auto lower_bound = curr_root->down;
        while (lower_bound != nullptr) {
            if (lower_bound->has_mid()) {
                double midpoint = lower_bound->get_mid();
                if (midpoint >= rect.lowY()) {
                    rect_total += lower_bound->up == nullptr ? 0.0 : lower_bound->up->measure_interval(rect.upX(), rect.lowX(), a, b);
                    lower_bound = lower_bound->down;
                } else {
                    lower_bound = lower_bound->up;
                }
            } else {
                rect_total += lower_bound->bottom_y >= rect.lowY() ? lower_bound->measure_interval(rect.upX(), rect.lowX(), a, b) : 0.0;
                break;
            }
        }
        return rect_total;
    }


    std::tuple<epoint_list_t, epoint_list_t> SlabTree::dyadic_approx(size_t upper, size_t lower) const {
        epoint_list_t mpts;
        epoint_list_t bpts;
        auto curr_root = get_containing(upper, lower);
        if (curr_root == nullptr) {
            return {mpts, bpts};
        }
        auto upper_bound = curr_root->up;
        while (upper_bound != nullptr) {
            if (upper_bound->has_mid()) {
                double midpoint = upper_bound->get_mid();
                if (midpoint < upper) {
                    if (upper_bound->down != nullptr) {
                        std::copy(upper_bound->down->m_merges.begin(), upper_bound->down->m_merges.end(), std::back_inserter(mpts));
                        std::copy(upper_bound->down->b_merges.begin(), upper_bound->down->b_merges.end(), std::back_inserter(bpts));
                    }
                    upper_bound = upper_bound->up;
                } else {
                    upper_bound = upper_bound->down;
                }
            } else {
                if (upper_bound->top_y <= upper) {
                    std::copy(upper_bound->m_merges.begin(), upper_bound->m_merges.end(), std::back_inserter(mpts));
                    std::copy(upper_bound->b_merges.begin(), upper_bound->b_merges.end(), std::back_inserter(bpts));
                }
                break;
            }
        }
        auto lower_bound = curr_root->down;
        while (lower_bound != nullptr) {
            if (lower_bound->has_mid()) {
                double midpoint = lower_bound->get_mid();
                if (midpoint >= lower) {
                    if (lower_bound->up != nullptr) {
                        std::copy(lower_bound->up->m_merges.begin(), lower_bound->up->m_merges.end(), std::back_inserter(mpts));
                        std::copy(lower_bound->up->b_merges.begin(), lower_bound->up->b_merges.end(), std::back_inserter(bpts));
                    }
                    lower_bound = lower_bound->down;
                } else {
                    lower_bound = lower_bound->up;
                }
            } else {
                if (lower_bound->bottom_y >= lower) {
                    std::copy(lower_bound->m_merges.begin(), lower_bound->m_merges.end(), std::back_inserter(mpts));
                    std::copy(lower_bound->b_merges.begin(), lower_bound->b_merges.end(), std::back_inserter(bpts));
                }
                break;
            }
        }
        std::sort(mpts.begin(), mpts.end());
        std::sort(bpts.begin(), bpts.end());
        return {mpts, bpts};
    }


    std::tuple<ERectangle, double> SlabTree::max_rectangle_slow(double m_a, double b_b) {

        ERectangle maximum_rectangle;
        double max_v = 0;
        std::vector<size_t> vert_splits = get_vert_breakpoints();
        for (size_t i = 0; i < vert_splits.size() - 1; i++) {
            for (size_t j = i + 1; j < vert_splits.size(); j++) {
                auto [mpts, bpts] = dyadic_approx(vert_splits[j], vert_splits[i]);

                std::vector<size_t> indices;
                std::vector<double> weights;
                for (auto& mpt : mpts) {
                    indices.emplace_back(mpt.get_x());
                    weights.emplace_back(mpt.get_weight() * m_a / this->total_m);
                }
                for (auto& bpt : bpts) {
                    indices.emplace_back(bpt.get_x());
                    weights.emplace_back(bpt.get_weight() * b_b / this->total_b);
                }

                auto interval = max_interval(indices, weights);
                if (interval.get_v() > max_v) {
                    maximum_rectangle = ERectangle(interval.get_r() + 1, vert_splits[j], interval.get_l(), vert_splits[i]);
                    max_v = interval.get_v();
                }
            }
        }
        return std::make_tuple(maximum_rectangle, max_v);
    }

    std::vector<MaxIntervalAlt> insert_updates(std::vector<MaxIntervalAlt> const& max_intervals,
                                              epoint_list_t const& updates, double scale) {
        std::vector<MaxIntervalAlt> new_set;
        for (auto& p : updates) {
            new_set.emplace_back(p.get_x(), p.get_weight() * scale);
        }
        std::vector<MaxIntervalAlt> updated;
        std::merge(max_intervals.begin(), max_intervals.end(), new_set.begin(), new_set.end(),
                std::back_inserter(updated), [](MaxIntervalAlt const& m1, MaxIntervalAlt const& m2) {
           return m1.get_l() < m2.get_l();
        });
        return updated;
    }

    std::vector<MaxIntervalAlt> update_mx_intervals(std::vector<MaxIntervalAlt> const & max_intervals,
            slab_ptr slab,
            double m_a, double b_b) {
        if (slab != nullptr) {
            return insert_updates(insert_updates(max_intervals, slab->m_merges, m_a), slab->b_merges, b_b);
        } else {
            return max_intervals;
        }
    }

    std::vector<MaxIntervalAlt> reduce_merges(std::vector<MaxIntervalAlt> const& max_intervals,
                                        std::vector<size_t> const& curr_splits) {

        if (max_intervals.empty()) {
            return {};
        }
        size_t j = 0;
        long prev = -1;
        std::vector<MaxIntervalAlt> merged;
        merged.reserve(max_intervals.size());
        merged.emplace_back(max_intervals.front());
        for (size_t i = 0; max_intervals.size() != 0 && i < max_intervals.size() - 1;) {
            if (prev < static_cast<long>(max_intervals[i].get_l()) && (curr_splits.size() <= j ||
                    max_intervals[i + 1].get_r() < curr_splits[j])) {
                merged.back() += max_intervals[i + 1];
                ++i;
            } else if (j < curr_splits.size() && curr_splits[j] <= max_intervals[i + 1].get_r()) {
                prev = static_cast<long>(curr_splits[j]);
                j++;
            } else {
                merged.emplace_back(max_intervals[i + 1]);
                i++;
            }
        }
        return merged;
    }


    slab_ptr get_top(slab_ptr p) {
        return p == nullptr ? nullptr : p->up;
    }
    slab_ptr get_bottom(slab_ptr p) {
        return p == nullptr ? nullptr : p->down;
    }



    using slab_frame = std::tuple<slab_ptr, slab_ptr, slab_ptr, slab_ptr, size_t,
                                 size_t,
                                 std::vector<MaxIntervalAlt>>;

    void emplace_helper(std::vector<slab_frame>& slabs, slab_ptr up, slab_ptr down, size_t high, size_t low,
            std::vector<MaxIntervalAlt> const& ms) {
        if (up != nullptr) {
            high = up->bottom_y;
        }
        if (down != nullptr) {
            low = down->top_y;
        }
        slabs.emplace_back(get_top(up), get_bottom(up),
                get_top(down), get_bottom(down),
                high, low, ms);
    }

    std::vector<MaxIntervalAlt> reduce_merges_helper(std::vector<MaxIntervalAlt> const& max_intervals,
                                              slab_ptr top, slab_ptr bottom) {

        if (top != nullptr  && bottom != nullptr) {
            std::vector<size_t> splits;
            std::set_union(top->global_split_offset.begin(), top->global_split_offset.end(),
                       bottom->global_split_offset.begin(), bottom->global_split_offset.end(),
                       std::back_inserter(splits));
            auto tmp_merges = reduce_merges(max_intervals, splits);
            return tmp_merges;
        } else if (top != nullptr) {
            return reduce_merges(max_intervals, top->global_split_offset);
        } else if (bottom != nullptr) {
            return reduce_merges(max_intervals, bottom->global_split_offset);
        } else {
            std::vector<size_t> splits;
            return reduce_merges(max_intervals, splits);
        }
    }


    std::tuple<ERectangle, double> SlabTree::max_rectangle_midpoint(double m_a, double b_b) {
        //Initialize with list of maximum intervals.
        ERectangle max_rect;
        m_a = m_a / total_m;
        b_b = b_b / total_b;

        std::vector<slab_frame> slab_stack;
        slab_stack.emplace_back(get_top(root->up), get_bottom(root->up),
                get_top(root->down), get_bottom(root->down),
                root->get_mid(), root->get_mid(), std::vector<MaxIntervalAlt>());

        double max_v = 0;
        while (!slab_stack.empty()) {
            auto [top_top, top_bottom, bottom_top, bottom_bottom, p_up, p_low, m4] = slab_stack.back();
            slab_stack.pop_back();
            if (!(top_top == nullptr &&
                  top_bottom == nullptr &&
                  bottom_top == nullptr &&
                  bottom_bottom == nullptr)) {
                // 4 doesn't have any merges
                // 3 is top_bottom merged.
                // 2 is bottom_top merged.
                // 1 is top_bottom and bottom_top merged.
                //Pretty sure this works.
                auto m3 = update_mx_intervals(m4, top_bottom, m_a, b_b);
                auto m1 = update_mx_intervals(m3, bottom_top, m_a, b_b);
                auto m2 = update_mx_intervals(m4, bottom_top, m_a, b_b);

                //Pretty sure this works.
                m1 = reduce_merges_helper(m1, top_top, bottom_bottom);
                m2 = reduce_merges_helper(m2, top_bottom, bottom_bottom);
                m3 = reduce_merges_helper(m3, top_top, bottom_top);
                m4 = reduce_merges_helper(m4, top_bottom, bottom_top);

                emplace_helper(slab_stack, top_top, bottom_bottom, p_up, p_low, m1);
                emplace_helper(slab_stack, top_bottom, bottom_bottom, p_up, p_low, m2);
                emplace_helper(slab_stack, top_top, bottom_top, p_up, p_low, m3);
                emplace_helper(slab_stack, top_bottom, bottom_top, p_up, p_low, m4);

            } else {
                if (!m4.empty() && max_v < m4[0].get_max().get_v()) {
                    size_t lx = m4[0].get_max().get_l();
                    size_t rx = m4[0].get_max().get_r();

                    max_rect = ERectangle(rx + 1, p_up, lx, p_low);
                    max_v = m4[0].get_max().get_v();
                }
            }
        }
        return std::make_tuple(max_rect, max_v);

    }

    std::tuple<ERectangle, double> SlabTree::max_rectangle(double m_a, double b_b) {
        std::vector<SlabTree> child_instances;
        child_instances.emplace_back(*this);
        ERectangle max_rect;
        double max_v = 0;
        while (!child_instances.empty()) {
            auto curr_tree = child_instances.back();
            child_instances.pop_back();
            auto [erect, eval] = curr_tree.max_rectangle_midpoint(m_a, b_b);
            if (eval > max_v) {
                max_rect = erect;
                max_v = eval;
            }
            if (curr_tree.get_root()->up != nullptr) {
                child_instances.emplace_back(curr_tree.get_upper_tree());
            }
            if (curr_tree.get_root()->down != nullptr) {
                child_instances.emplace_back(curr_tree.get_lower_tree());
            }
        }
        return std::make_tuple(max_rect, max_v);
    }


    std::tuple<Rectangle, double> max_rectangle(const wpoint_list_t& mpts, const wpoint_list_t& bpts, double eps, double a, double b) {
        auto [m_pts, b_pts, xmap, ymap] = pyscan::to_epoints(mpts, bpts);
        SlabTree tree(m_pts, b_pts, eps);
        tree.even_compress(eps / log(1 / eps));
        tree.compute_splits();
        auto [max_rect, max_v] = tree.max_rectangle(a, b);
        return std::make_tuple(Rectangle(xmap[max_rect.upX()], ymap[max_rect.upY()], xmap[max_rect.lowX()], ymap[max_rect.lowY()]), max_v);
    }

//    std::tuple<Rectangle, double> max_rectangle_convex(wpoint_list_t const& mpts, wpoint_list_t const& bpts, double eps, discrepancy_func_t const &f) {
//        /*
//         * This uses the
//         */
//        auto [m_pts, b_pts, xmap, ymap] = pyscan::to_epoints(mpts, bpts);
//        auto phi = [&] (Vec2 const& v) {return f(v[0], 1.0,  v[1], 1.0); };
//        SlabTree tree(m_pts, b_pts, 1 / eps);
//
//        double maxV = 0;
//        auto linemaxF = [&] (Vec2 const& dir) {
//            auto [erect, val] = tree.max_rectangle(dir[0], dir[1]);
//
//
//            double curr_val = phi( Vec2{val * dir[0], val * dir[1] });
//            if (curr_val > maxV) {
//                max_subgrid = curr_subgrid;
//                maxV = curr_val;
//            }
//            return curr_mb;
//        };
//        approximateHull(eps, phi, linemaxF);
//        return max_subgrid;
//    }
}
