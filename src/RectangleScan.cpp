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
            tmp(0) = i;

            x_map[i] = access(mpts, bpts, ixx[i])(0);
            y_map[i] = access(mpts, bpts, ixy[i])(1);

            tmp.set_weight(access(mpts, bpts, ixx[i]).get_weight());
            access(empts, ebpts, ixy[i])(1) = i;
        }
        return std::make_tuple(empts, ebpts, x_map, y_map);
    }

     MaxIntervalAlt::MaxIntervalAlt(size_t val, double weight) : left_max(val, val, weight), right_max(val, val, weight), center_max(val, val, weight) {}
     MaxIntervalAlt::MaxIntervalAlt(size_t lpt, size_t rpt) : left_max(lpt, rpt, 0.0), right_max(lpt, rpt, 0.0), center_max(lpt, rpt, 0.0) {}

     MaxIntervalAlt& MaxIntervalAlt::operator+=(const MaxIntervalAlt& op) {
         if (center_max.get_v() < op.center_max.get_v()) {
             center_max = op.center_max;
         }

         if (right_max.get_v() + op.left_max.get_v() > center_max.get_v()) {
             center_max = Interval(right_max.get_l(), op.left_max.get_r(), right_max.get_v() + op.left_max.get_v());
         }

         //If left max is whole interval then we might need to extend the interval.
         if (left_max.get_r() == right_max.get_r() && left_max.get_v() + op.left_max.get_v() > left_max.get_v()) {
             left_max = Interval(left_max.get_l(), op.left_max.get_r(), left_max.get_v() + op.left_max.get_v());
         }

         // op.right max is whole interval of op then we might need to extend op.right_max to include right_max
         if (op.right_max.get_l() == op.left_max.get_l() && op.right_max.get_v() + right_max.get_v() > op.right_max.get_v()) {
             right_max = Interval(right_max.get_l(), op.right_max.get_r(), right_max.get_v() + op.right_max.get_v());
         } else {
             right_max = op.right_max;
         }
         return *this;
     }

     double MaxIntervalAlt::left() const {
         return left_max.get_l();
     }
     double MaxIntervalAlt::right() const {
         return right_max.get_r();
     }

     Interval MaxIntervalAlt::get_max() const {
         return center_max;
     }



     epoint_list_t compress(epoint_it_t b, epoint_it_t e, bool compression, double m_w) {
        /*
         * Takes a vector of weighted points and constructed a smaller set of weighted points equally spaced with
         * the a different set of weights.
         */
        if (compression) {
            double curr_w = 0;
            std::vector<ept_t> break_pts;
            for (; b != e; ++b) {
                //std::cout << m_w << " "  << curr_w << std::endl;

                if ((curr_w + b->get_weight()) > m_w) {
                    break_pts.emplace_back(*b);
                    curr_w = 0;
                } else {
                    curr_w += b->get_weight();
                }
            }
            return break_pts;
        } else {
             epoint_list_t merges;
             for (; b != e; ++b) {
                 merges.emplace_back(*b);
             }
             return merges;
         }
    }

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
        double curr_weight = 0;
        for (auto &p : merge_buffer) {
            curr_weight += p.get_weight();
            if (curr_weight >= max_w * (total_m + total_b)) {
                curr_weight = 0;
                vert_decomp.emplace_back(p.get_y());
            }
        }

        SlabTree::init(mpts, bpts, vert_decomp, true, max_w);
    }

    SlabTree::SlabTree(std::vector<size_t> const &vert_decomp, epoint_list_t mpts, epoint_list_t bpts, bool compression, double max_w) :
        total_m(computeTotal(mpts)), total_b(computeTotal(bpts)) {
        /*
         * Create a vertical decomposition of the point set so that we can split the points into a sequence of horizontal
         * strips where each contains at most max_w.
         */
        SlabTree::init(std::move(mpts), std::move(bpts), vert_decomp, compression, max_w);
    }

    void SlabTree::init(epoint_list_t mpts, epoint_list_t bpts, std::vector<size_t> const &vert_decomp, bool compression, double max_w) {

        //std::sort(vert_decomp.begin(), vert_decomp.end());

        using wrng_it = std::tuple<epoint_it_t , epoint_it_t >;
        using vert_it = std::vector<size_t>::const_iterator;
        using vrng_it = std::tuple<vert_it, vert_it>;

        using cell_t = std::tuple<slab_ptr, slab_ptr *, vrng_it, wrng_it, wrng_it>;
        using cell_list_t = std::vector<cell_t>;

        {
            double v_b = *(vert_decomp.begin());
            double v_e = *(vert_decomp.end() - 1);
            //Needs to match ERectangle->contains
            auto in_interval = [v_b, v_e] (const ept_t& pt) {
                return v_b <= pt(1) && pt(1) < v_e;
            };
            auto m_new_end = std::partition(mpts.begin(), mpts.end(), in_interval);
            mpts.erase(m_new_end, mpts.end());
            auto b_new_end = std::partition(bpts.begin(), bpts.end(), in_interval);
            bpts.erase(b_new_end, bpts.end());
        }

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
            //std::cout << m_e - m_b<< std::endl;
            auto m_merges = compress(m_b, m_e, compression, max_w * total_m);
            auto b_merges = compress(b_b, b_e, compression, max_w * total_b);

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
            } else {
                roots.emplace(*el_ptr);
            }
        }

        using s_it = std::vector<size_t>::iterator;
        auto inplace_set_union = [](std::vector<size_t>& memory, s_it b, s_it mid, s_it e) {
            // Now merge into the already existing split offsets.
            std::inplace_merge(b, mid, e);
            auto end_it = std::unique(b, e);
            memory.erase(end_it, e);
        };
        //Now go back through and create a list of all the lower splits.
        while (!roots.empty()) {
            auto child = roots.front();
            roots.pop();
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
                    roots.emplace(p);
                }
            }
        }
    }

    slab_ptr SlabTree::get_containing(ERectangle const &rect) const {
        assert(rect.upY() > rect.lowY());
        auto curr_root = root;
        while (curr_root != nullptr) {
            //std::cout << "Cur Root "<<  curr_root->top_y << " " << rect.upY() <<  " " << rect.lowY() << " " << curr_root->bottom_y << std::endl;
            //Check if this region is contained completely in the left or right branch.
            if (curr_root->has_mid()) {
                size_t mid = curr_root->get_mid();

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

    double SlabTree::measure_rect(ERectangle const &rect, double a, double b) const {
        auto curr_root = get_containing(rect);
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


    std::vector<MaxIntervalAlt> insert_updates(std::vector<MaxIntervalAlt> const& max_intervals,
                                              epoint_list_t const& updates, double scale) {
        std::vector<MaxIntervalAlt> new_set;
        for (auto& p : updates) {
            new_set.emplace_back(p.get_x(), p.get_weight() * scale);
        }
        std::vector<MaxIntervalAlt> updated;
        std::merge(max_intervals.begin(), max_intervals.end(), new_set.begin(), new_set.end(),
                std::back_inserter(updated), [](MaxIntervalAlt const& m1, MaxIntervalAlt const& m2) {
           return m1.left() < m2.left();
        });
        for (size_t i = 0; !updated.empty() && i < updated.size() - 1; i++){
            assert(updated[i].right() < updated[i + 1].left());
        }
        return updated;
    }

    std::vector<MaxIntervalAlt> update_mx_intervals(std::vector<MaxIntervalAlt> const & max_intervals, slab_ptr slab,
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
        merged.emplace_back(max_intervals.front());
        for (size_t i = 0; max_intervals.size() != 0 && i < max_intervals.size() - 1;) {
            if (prev < static_cast<long>(max_intervals[i].left()) && (curr_splits.size() <= j || max_intervals[i + 1].right() < curr_splits[j])) {
                merged.back() += max_intervals[i + 1];
                ++i;
            } else if (j < curr_splits.size() && curr_splits[j] <= max_intervals[i + 1].right()) {
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
            high = up->top_y;
        }
        if (down != nullptr) {
            low = down->bottom_y;
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
//            std::cout << "Before " << max_intervals << std::endl;
//            std::cout << "Splits " << splits << std::endl;
            auto tmp_merges = reduce_merges(max_intervals, splits);
//            std::cout << "After " << tmp_merges << std::endl;
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

    std::tuple<ERectangle, double> SlabTree::max_rectangle(double m_a, double b_b) {
        //Initialize with list of maximum intervals.
        ERectangle max_rect;
        m_a = m_a / total_m;
        b_b = b_b / total_b;

        std::vector<slab_frame> slab_stack;
        slab_stack.emplace_back(get_top(root->up), get_bottom(root->up),
                get_top(root->down), get_bottom(root->down),
                root->top_y, root->bottom_y, std::vector<MaxIntervalAlt>());

        double max_v = 0;
        while (!slab_stack.empty()) {
            auto [top_top, top_bottom, bottom_top, bottom_bottom, p_up, p_low, m4] = slab_stack.back();
//            std::cout << m4 << std::endl;
//            //std::cout << p_splts << std::endl;
//            std::cout << "[" <<  p_up << ", " << p_low << "] " << std::endl;
//            //std::cout << max_intervals << std::endl;
//            //std::cout << "Slab" << top_top << " " << top_bottom << " " << bottom_top << " " << bottom_bottom << std::endl;
//            std::cout << std::endl;
            slab_stack.pop_back();
            if (!(top_top == nullptr &&
                  top_bottom == nullptr &&
                  bottom_top == nullptr &&
                  bottom_bottom == nullptr)) {
                // 4 doesn't have any merges
                // 3 is top_bottom merged.
                // 2 is bottom_top merged.
                // 1 is top_bottom and bottom_top merged.
                auto m3 = update_mx_intervals(m4, top_bottom, m_a, b_b);
                auto m1 = update_mx_intervals(m3, bottom_top, m_a, b_b);
                auto m2 = update_mx_intervals(m4, bottom_top, m_a, b_b);

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


    std::tuple<Rectangle, double> max_rectangle(const wpoint_list_t& mpts, const wpoint_list_t& bpts, double eps, double a, double b) {
        auto [m_pts, b_pts, xmap, ymap] = pyscan::to_epoints(mpts, bpts);
        SlabTree tree(m_pts, b_pts, eps);
        //std::cout << tree << std::endl;
        std::vector<SlabTree> child_instances;
        child_instances.emplace_back(tree);
        ERectangle max_rect;
        double max_v = 0;
        while (!child_instances.empty()) {
            auto curr_tree = child_instances.back();
            child_instances.pop_back();
            auto [erect, eval] = curr_tree.max_rectangle(a, b);
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
