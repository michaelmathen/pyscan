//
// Created by mmath on 5/28/17.
//
#include <functional>
#include <tuple>
#include <unordered_set>
#include <zlib.h>

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


    Subgrid max_subgrid_convex_theory(Grid const &grid, double eps,
                                      discrepancy_func_t const &f) {
        /*
         * This uses the
         */
        auto phi = [&] (Vec2 const& v) {return f(v[0], 1.0, v[1], 1.0); };
        Subgrid curr_subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
        Subgrid max_subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
        double maxV = 0;
        auto linemaxF = [&] (Vec2 const& dir) {
            curr_subgrid = max_subgrid_linear_theory(grid, grid.size() * r_ratio(grid.size()), dir[0], dir[1]);
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

//                //Remove duplicates in each cell.
//                // this isn't necessary, but will hopefully speed things up.
//                for (auto& cell : cells) {
//                    std::sort(cell.begin(), cell.end(), [](LabelW const& l1, LabelW const& l2) {
//                        return l1.label < l2.label;
//                    });
//                    auto it = std::unique(cell.begin(), cell.end(), []( LabelW const& l1, LabelW const& l2) {
//                        return l1.label == l2.label;
//                    });
//                    cell.resize(static_cast<unsigned long>(std::distance(cell.begin(), it)));
//
//                }
            };
            grid_insert(m_labels, m_points);
            grid_insert(b_labels, b_points);
//            std::cout << x_values << std::endl;
//            std::cout << y_values << std::endl;
//
//            std::cout << x_values.size() << std::endl;
//            std::cout << y_values.size() << std::endl;
//
//            std::cout << m_labels.size() << std::endl;
//            std::cout << b_labels.size() << std::endl;
//
//            std::cout << m_labels << std::endl;
//            std::cout << b_labels << std::endl;

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



    ////////////////////////////////////////////////////////////////////////////////
    //Maximum rectangle from paper code/////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    ValueInterval::ValueInterval() : left(0), value(0), right(0) {}
    ValueInterval::ValueInterval(double val, size_t left, size_t right) : left(left),
                                                               value(val),
                                                               right(right)
   {}
    void ValueInterval::print(std::ostream& os) const {
        os << "[" << left << ", " << right << "]";
        os << " = " << value;
    }
    size_t ValueInterval::getLeft() const { return left; }
    size_t ValueInterval::getRight() const { return right; }
    void ValueInterval::setLeft(size_t left_c) { left = left_c; }
    void ValueInterval::setRight(size_t right_c) { right = right_c; }
    double ValueInterval::getValue() const {
        return value;
    }

    void ValueInterval::setValue(double val) {
        value = val;
    }

    ValueInterval &ValueInterval::operator+=(double val) {
        value += val;
        return *this;
    }


    MaxInterval::MaxInterval() : left_max(), right_max(), center_max() {}
    MaxInterval::MaxInterval(double value, size_t index) : left_max(value, index, index),
                                              right_max(value, index, index),
                                              center_max(value, index, index) {}


   MaxInterval &MaxInterval::operator+=(MaxInterval const &op) {

        ValueInterval merged_max;
        //Figure out the new center max value
        if (right_max.getValue() < 0 || op.left_max.getValue() < 0) {
            if (right_max.getValue() < op.left_max.getValue()) {
                merged_max = op.left_max;
            } else {
                merged_max = right_max;
            }
        } else {
            merged_max = ValueInterval(right_max.getValue() + op.left_max.getValue(),
                                       right_max.getLeft(), op.left_max.getRight());
        }

        if (merged_max.getValue() >= center_max.getValue() && merged_max.getValue() >= op.center_max.getValue()) {
            center_max = merged_max;
        } else if (op.center_max.getValue() > center_max.getValue()) {
            center_max = op.center_max;
        } // else it remains the same.
        //Figure out the new left max value
        if (left_max.getRight() == this->getRight() && op.left_max.getValue() >= 0) {
            left_max.setRight(op.left_max.getRight());
            left_max.setValue(op.left_max.getValue() + left_max.getValue());
        }
        //Figure out the new right max value
        if (op.right_max.getLeft() == op.getLeft() && right_max.getValue() >= 0) {
            //Left boundary of the right max doesn't change
            right_max.setValue(op.right_max.getValue() + right_max.getValue());
        } else {
            // New left boundary of the right max.
            right_max = op.right_max;
        }
        right_max.setRight(op.getRight());
        return *this;
    }

    size_t MaxInterval::getLeft() const { return left_max.getLeft(); }
    size_t MaxInterval::getRight() const { return right_max.getRight(); }
    size_t MaxInterval::lowCIx() const {
        return center_max.getLeft();
    }

    size_t MaxInterval::upCIx() const {
        return center_max.getRight();
    }

    double MaxInterval::lValue() const {
        return left_max.getValue();
    }

    double MaxInterval::rValue() const {
        return right_max.getValue();
    }

    double MaxInterval::cValue() const {
        return center_max.getValue();
    }

    void MaxInterval::print(std::ostream& os) const {
        os << "[" << getLeft() << ", " << getRight() << "]";
        os << "= {";
        os << left_max << ", " << center_max << ", " << right_max << "}";
    }


    MaximumIntervals::MaximumIntervals(size_t l, size_t r) : left_col(l), right_col(r) {}

    void MaximumIntervals::print(std::ostream& os) const {
        os << "[" << left_col << ", " << right_col << "]";
        os << max_intervals << std::endl;
        os << weights << std::endl;
        os << weight_index << std::endl;

        //os << intervals << std::endl;
    }

    size_t MaximumIntervals::getWeightNum() const {
        return weights.size();
    }

    size_t MaximumIntervals::getIntervalNum() const {
        return intervals.size();
    }

    MaximumIntervals::MaximumIntervals(size_t r) : intervals(r, I_Type::VALUE),
                                 max_intervals(),
                                 weights(r, 0),
                                 weight_index(r, 0),
                                 left_col(0),
                                 right_col(r - 1) {
        for (size_t i = 0; i < r; i++)
            weight_index[i] = i;
    }

    void MaximumIntervals::setBounds(size_t l_c, size_t r_c) {
        left_col = l_c;
        right_col = r_c;
    }

    void MaximumIntervals::updateBounds(size_t l_c, size_t r_c) {
        left_col = std::min(l_c, left_col);
        right_col = std::max(r_c, right_col);
    }

    /*
     * Indices and weights are assumed to be presorted.
     */
    void MaximumIntervals::updateWeights(std::vector<double> const &new_weights,
                       std::vector<size_t> const &indices,
                       double w) {

        size_t j = 0;
        for (size_t i = 0; i < new_weights.size(); i++) {
            while (j < weight_index.size() && weight_index[j] <= indices[i]) {
                if (weight_index[j] == indices[i]) {
                    weights[j] += new_weights[i] * w;
                }
                j++;
            }
        }
    }

    MaximumIntervals MaximumIntervals::mergeZeros(BloomFilter const& f1) const {
        auto mightBePresent = [&f1](uint32_t index) {
            return f1.mightBePresent(index);
        };

        return mergeZeros(mightBePresent);
    }
    MaximumIntervals MaximumIntervals::mergeZeros(BloomFilter const& f1, BloomFilter const& f2) const {
        auto mightBePresent = [&f1, &f2](uint32_t index) {
            return f1.mightBePresent(index) || f2.mightBePresent(index);
        };
        return mergeZeros(mightBePresent);
    }

    Subgrid MaximumIntervals::getMax() const {
        if (!intervals.empty()) {
            int mx_ix = 0;
            int w_ix = 0;
            MaxInterval new_max;
            if (intervals[0] == I_Type::VALUE) {
                new_max = MaxInterval(weights[0], weight_index[0]);
                w_ix += 1;
            } else {
                new_max = max_intervals[0];
                mx_ix += 1;
            }
            for (size_t i = 1; i < intervals.size(); i++) {
                if (intervals[i] == I_Type::VALUE) {
                    new_max += MaxInterval(weights[w_ix], weight_index[w_ix]); // Merge the previous
                    w_ix += 1;
                } else {
                    new_max += max_intervals[mx_ix]; // Merge the previous
                    mx_ix += 1;
                }
            }
            return {new_max.upCIx(), this->getRight(), new_max.lowCIx(), this->getLeft(), new_max.cValue()};
        } else {
            return {this->getLeft(), 0, this->getRight(), 0, -std::numeric_limits<double>::infinity()};
        }
    }

    size_t MaximumIntervals::maxIntervalNum() const {
        return max_intervals.size();
    }

    size_t MaximumIntervals::valueNum() const {
        return weights.size();
    }

    MaxInterval MaximumIntervals::getFirstMax() const {
        return max_intervals[0];
    }

    size_t MaximumIntervals::getLeft() const {
        return left_col;
    }

    size_t MaximumIntervals::getRight() const {
        return right_col;
    }

    void MaximumIntervals::setLeft(size_t left) {
        left_col = left;
    }

    void MaximumIntervals::setRight(size_t right) {
        right_col = right;
    }


    SlabNode::SlabNode(Grid const& g, size_t left_col, size_t right_col, size_t r_prime) :
            left_col(left_col),
            right_col(right_col) {
        double red_weight = 0;
        double blue_weight = 0;
        //std::cout << "[" << left_col << ", " << right_col << "]" << std::endl;
        for (size_t row = 0; row < g.size(); row++) {
            for (size_t col = left_col; col <= right_col; col++) {
                red_weight += g.redWeight(row, col);
                blue_weight += g.blueWeight(row, col);
            }
            if (red_weight + blue_weight > 1 / static_cast<double>(r_prime)) {
                red_weights.push_back(red_weight);
                blue_weights.push_back(blue_weight);
                red_weight = 0;
                blue_weight = 0;
                indices.push_back(row);
            }
        }
        //std::cout << red_weights << std::endl;
        //std::cout << indices << std::endl;
    }

    size_t SlabNode::getMid() const {
        return (left_col + right_col) / 2;
    }

    /*
    double measureInterval(size_t l_r, size_t u_r, double a, double b) const {
       double sum = 0;
       //std::cout << l_r << " " << u_r << std::endl;
       //std::cout << indices << std::endl;
       for (size_t i = 0; i < indices.size(); i++) {
           if (l_r <= indices[i] && indices[i] <= u_r) {
               sum += a * red_weights[i] + b * blue_weights[i];
           }
        }
        return sum;
    }


    double measureRect(size_t l_c, size_t l_r, size_t u_c, size_t u_r, double a, double b) const {
        // If it completely overlaps
        auto measureRectH = [&](slab_ptr ptr) {
            if (ptr != nullptr) {
                return ptr->measureRect(l_c, l_r, u_c, u_r, a, b);
            } else {
                return 0.0;
            }
        };
        //std::cout << l_c << " " << l_r << " " << u_c << " " << u_r << std::endl;
        //std::cout << left_col << " " << right_col << std::endl;
        //print(std::cout);
        //std::cout << std::endl;

        if (l_c <= left_col && right_col <= u_c) {
            return measureInterval(l_r, u_r, a, b);
        } else if (right_col < l_c || u_c < left_col) {
            return 0.0;
        } else {
            return measureRectH(left) + measureRectH(right);
        }
    }
    */
    bool SlabNode::hasChildren() const {
        return left != nullptr;
    }

    void SlabNode::print(std::ostream& os) const {
        os << "red= " << red_weights << std::endl;
        os << "blue= " << blue_weights << std::endl;
        os << "indices= " << indices << std::endl;
    }


    using slab_ptr = SlabNode::slab_ptr;

    struct SlabTree {
        size_t r;
        slab_ptr root;
        /*
         * red_b, red_e -- iterators to the red points
         * blue_b, blue_e --iterators to the blue points
         * r -- the number of points to use in the top level slab.
         */
        SlabTree(Grid const &grid, size_t r_prime) : r(grid.size()), root(nullptr) {
            size_t upper = r - 1, lower = 0;
            std::deque<std::tuple<int, int, slab_ptr>> stack_nodes;
            root = std::make_shared<SlabNode>(grid, lower, upper, r_prime);
            stack_nodes.emplace_back(lower, upper, root);
            while (!stack_nodes.empty()) {
                auto node = stack_nodes.back();
                if (std::get<1>(node) <= std::get<0>(node)) {
                    stack_nodes.pop_back();
                    continue;
                }

                auto mid = (std::get<1>(node) + std::get<0>(node)) / 2;

                auto &parent = std::get<2>(node);
                auto left_child = std::make_shared<SlabNode>(grid, std::get<0>(node), mid, r_prime);
                auto right_child = std::make_shared<SlabNode>(grid, mid + 1, std::get<1>(node), r_prime);
                parent->left = left_child;
                parent->right = right_child;
                stack_nodes.pop_back();
                stack_nodes.emplace_back(std::get<0>(node), mid, left_child);
                stack_nodes.emplace_back(mid + 1, std::get<1>(node), right_child);
            }

            computeBlooms(root);
        }

        void print(std::ostream& os) const {
            printHelper(root, os);
        }
        /*
        double measure(size_t l_c, size_t l_r, size_t u_c, size_t u_r, double a, double b) const {
            if (root == nullptr) {
                return 0.0;
            } else {
                return root->measureRect(l_c, l_r, u_c, u_r, a, b);
            }
        }

        double measureSubgrid(Subgrid const& subgrid, double a, double b) const {
            return measure(subgrid.upX(), subgrid.lowRow(), subgrid.upCol(), subgrid.upRow(), a, b);
        }
         */

    private:

        std::vector<int> computeBlooms(slab_ptr &curr_root) {
#ifdef DEBUG
            std::cout << "Computing blooms in slab tree" << std::endl;
#endif
            if (curr_root == nullptr) {
                return std::vector<int>();
            } else {
                //curr_root->print(std::cout);
                std::vector<int> tmp_indices;
                std::vector<int> new_indices;
                auto left_non_zeros = computeBlooms(curr_root->left);
                auto right_non_zeros = computeBlooms(curr_root->right);
                //std::cout << "left = " << left_non_zeros << std::endl;
                // std::cout << "right = " << right_non_zeros << std::endl;
                auto &curr_indices = curr_root->indices;
                //std::cout << "current = " << curr_root->indices << std::endl;
                std::set_union(left_non_zeros.begin(), left_non_zeros.end(),
                               right_non_zeros.begin(), right_non_zeros.end(),
                               std::back_inserter(tmp_indices)
                );
                std::set_union(tmp_indices.begin(), tmp_indices.end(),
                               curr_indices.begin(), curr_indices.end(),
                               std::back_inserter(new_indices));
                //std::cout << "merged = " << new_indices << std::endl;
                curr_root->non_zeros = BloomFilter(new_indices.size(), .1);

                for (size_t i = 0; i < new_indices.size(); i++) {
                    //std::cout << curr_root->non_zeros.mightBePresent(i) << std::endl;
                    curr_root->non_zeros.insert(new_indices[i]);
                }
                return new_indices;
            }
        }

        void printHelper(slab_ptr curr_root, std::ostream& os) const {
            if (curr_root != nullptr) {
                curr_root->print(os);
                printHelper(curr_root->left, os);
                printHelper(curr_root->right, os);
            }
        }
    };

     void mergeWeights(MaximumIntervals &interval, slab_ptr const& left, slab_ptr const& right) {
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

    Subgrid maxSubgridSpan(slab_ptr left_side, slab_ptr right_side, MaximumIntervals const &maxIntInitial,
                             double a, double b) {

        using SubProblem = std::tuple<slab_ptr, slab_ptr, MaximumIntervals>;
        std::deque<SubProblem> sub_problem_stack;
        sub_problem_stack.emplace_back(left_side, right_side, maxIntInitial);
        Subgrid max_subgrid = Subgrid(-1, -1, -1, -1, -std::numeric_limits<double>::infinity());

        while (!sub_problem_stack.empty()) {
           SubProblem sub_problem = sub_problem_stack.back();
           sub_problem_stack.pop_back();

           auto& [left, right, maxInt] = sub_problem;
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


    Subgrid max_subgrid_lin(slab_ptr slab, MaximumIntervals const &maxInt, double a, double b) {
        if (slab == nullptr) {
            //std::cout << maxInt << maxInt.getLeft() << maxInt.getRight() <<  std::endl;
            return maxInt.getMax();
        } else {
            Subgrid max = Subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
            std::deque<std::tuple<slab_ptr, MaximumIntervals>> slabs;
            slabs.emplace_back(slab, maxInt);
            while (!slabs.empty()) {
                auto slab_mx = slabs.back();
                auto& [new_slab, mx] = slab_mx;
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
                    slabs.emplace_back(new_slab->left, new_maxInt);
                    slabs.emplace_back(new_slab->right, new_maxInt);
                }
            }
            return max;
        }
    }

    Subgrid max_subgrid_linear_theory(Grid const &grid, long r_prime, double a, double b) {
        SlabTree slabT(grid, (size_t)r_prime);
#ifdef DEBUG
        std::cout<<"finished the slab tree" << std::endl;
#endif
        //std::cout << slabT << std::endl;
        MaximumIntervals maxInt(grid.size());
        return max_subgrid_lin(slabT.root, maxInt, a, b);
    }

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


    std::tuple<std::vector<pt2_t*>, std::vector<double>> compress(const std::vector<wpt2_t*> & wpts, double m_w) {
        /*
         * Takes a vector of weighted points and constructed a smaller set of weighted points equally spaced with
         * the a different set of weights.
         */
        if (wpts.empty()) {
            return std::make_tuple(std::vector<pt2_t*>(), std::vector<double>());
        }
        double curr_weight = (*wpts.begin())->get_weight();
        std::vector<double> weights;
        std::vector<pt2_t *> break_pts;
        for (auto it = wpts.begin() + 1; it != wpts.end(); ++it) {
            if ((curr_weight + (*it)->get_weight()) > m_w) {
                weights.emplace_back(curr_weight);
                break_pts.emplace_back(*(it - 1));
                curr_weight = 0;
            }
            curr_weight += (*it)->get_weight();
        }
        if (curr_weight == m_w) {
            weights.emplace_back(curr_weight);
            break_pts.emplace_back(*(wpts.end() - 1));
        }
        return std::make_tuple(break_pts, weights);
    }

    void normalize(wpoint_list_t& pts) {
        //Compute the total weight.
        double total_w = 0;
        for (auto& p : pts)  total_w += p.get_weight();
        std::transform(pts.begin(), pts.end(), pts.begin(), [total_w](const wpt2_t& pt){
            return wpt2_t(pt.get_weight() / total_w, pt[0], pt[1], pt[2]);
        });
    };

    class SlabTreeAlt {

        class Slab {

        public:
            std::vector<pt2_t*> split_offsets;
            std::vector<pt2_t*> m_merges;
            std::vector<double> m_weights;
            std::vector<pt2_t*> b_merges;
            std::vector<double> b_weights;

            pt2_t* top_y;
            pt2_t* bottom_y;

            Slab* down = nullptr;
            Slab* up = nullptr;
            Slab* parent = nullptr;

            Slab(Slab* p,
                 std::vector<pt2_t*> m_m,
                 std::vector<double> m_w,
                 std::vector<pt2_t*> b_m,
                 std::vector<double> b_w) :
                                            m_merges(std::move(m_m)),
                                            m_weights(std::move(m_w)),
                                            b_merges(std::move(b_m)),
                                            b_weights(std::move(b_w)),
                                            parent(p) {}
        };
        std::vector<Slab> slab_heap;

        wpoint_list_t mpts;
        wpoint_list_t bpts;
    public:
        using slab_it_t = std::vector<Slab>::const_iterator;

        SlabTreeAlt(wpoint_list_t ms, wpoint_list_t bs, double max_w) : mpts(std::move(ms)), bpts(std::move(bs)) {
            normalize(mpts);
            normalize(bpts);

            using wptr_list_t = std::vector<wpt2_t*>;

            wptr_list_t net_pts;
            wptr_list_t m_sample_pts;
            m_sample_pts.reserve(mpts.size());
            for (auto& p : mpts) {
                m_sample_pts.emplace_back(&p);
                net_pts.emplace_back(&p);
            }
            wptr_list_t b_sample_pts;
            b_sample_pts.reserve(bpts.size());
            for (auto& p : bpts) {
                b_sample_pts.emplace_back(&p);
                net_pts.emplace_back(&p);
            }

            //Sort net verticaly
            std::sort(net_pts.begin(), net_pts.end(), [](wpt2_t* p1, wpt2_t* p2) {
                return p1->operator()(1) < p2->operator()(1);
            });
            //Sort samples horizontally
            auto x_order = [](pt2_t* p1, pt2_t* p2) {
                return p1->operator()(0) < p2->operator()(0);
            };
            std::sort(m_sample_pts.begin(), m_sample_pts.end(), x_order);
            std::sort(b_sample_pts.begin(), b_sample_pts.end(), x_order);

            std::vector<std::tuple<Slab*, wptr_list_t, wptr_list_t, wptr_list_t, double, bool> > active;
            active.emplace_back(nullptr, net_pts, m_sample_pts, b_sample_pts, 2.0, true);

            auto weighted_median_sorted = [max_w](wptr_list_t const& w_arr) {
                double curr_weight = 0;
                for (auto pt : w_arr) {
                    curr_weight += pt->get_weight();
                    if (curr_weight >= max_w) {
                        return pt;
                    }
                }
                return static_cast<wpt2_t *>(nullptr);
            };

            std::vector<Slab*> roots;
            while (!active.empty())  {
                auto [parent, local_net_pts, m_node_pts, b_node_pts, division, up] = active.back();
                active.pop_back();
                //Split m_pts and b_pts in this slab by the weighted median.
                auto split_pt = weighted_median_sorted(local_net_pts);

                auto net_iter_splt = std::partition(local_net_pts.begin(), local_net_pts.end(), [split_pt](wpt2_t* ptr){
                    return ptr->operator()(1) < split_pt->operator()(1);
                });

                auto m_iter_splt = std::partition(m_node_pts.begin(), m_node_pts.end(), [split_pt](wpt2_t* ptr){
                    return ptr->operator()(1) < split_pt->operator()(1);
                });

                auto b_iter_splt = std::partition(b_node_pts.begin(), b_node_pts.end(), [split_pt](wpt2_t* ptr){
                    return ptr->operator()(1) < split_pt->operator()(1);
                });

                //Compress the current m_pts and b_pts and compute the weights.
                auto [m_merges, m_weights] = compress(m_node_pts, max_w);
                auto [b_merges, b_weights] = compress(b_node_pts, max_w);
                //Compute the slab now.
                slab_heap.emplace_back(parent, m_merges, m_weights, b_merges, b_weights);
                //Update the parent with this child and set the correct boundaries.
                if (parent != nullptr) {
                    if (up) {
                        parent->up = &slab_heap.back();
                        slab_heap.back().top_y = parent->top_y;
                        slab_heap.back().bottom_y = split_pt;
                    } else {
                        parent->down = &slab_heap.back();
                        slab_heap.back().bottom_y = parent->bottom_y;
                        slab_heap.back().top_y = split_pt;
                    }
                } else {
                    slab_heap.back().top_y = local_net_pts.front();
                    slab_heap.back().bottom_y = local_net_pts.back();
                }
                //Now emplace back to the active
                if (1 / (2 * division) > max_w) {
                    active.emplace_back(&slab_heap.back(),
                            wptr_list_t(local_net_pts.begin(), net_iter_splt),
                            wptr_list_t(m_node_pts.begin(), m_iter_splt),
                            wptr_list_t(b_node_pts.begin(), b_iter_splt),
                            2 * division, false);
                    active.emplace_back(&slab_heap.back(),
                            wptr_list_t(net_iter_splt, local_net_pts.end()),
                            wptr_list_t(m_iter_splt, m_node_pts.end()),
                            wptr_list_t(b_iter_splt, b_node_pts.end()),
                            2 * division, true);
                } else {
                    roots.emplace_back(&slab_heap.back());
                }
            }

            //Now go back through and create a list of all the lower splits.
            while (!roots.empty()) {
                auto child = roots.back();
                roots.pop_back();

                //First merge the m_merges and b_merges
                child->split_offsets.reserve(child->split_offsets.size() + child->m_merges.size() + child->b_merges.size());

                auto mid = child->split_offsets.end();
                std::merge(child->m_merges.begin(),
                        child->m_merges.end(),
                        child->b_merges.begin(),
                        child->b_merges.end(),
                        std::back_inserter(child->split_offsets),
                        x_order);

                //Now merge into the already existing split offsets.
                std::inplace_merge(child->split_offsets.begin(), mid, child->split_offsets.end(), x_order);

                //Now merge into the parent.
                if (child->parent != nullptr) {
                    auto p = child->parent;
                    auto mid_p = p->split_offsets.end();
                    std::copy(child->split_offsets.begin(), child->split_offsets.end(),
                            std::back_inserter(p->split_offsets));
                    std::inplace_merge(p->split_offsets.begin(), mid_p, p->split_offsets.end(), x_order);
                }
            }

        }




        std::tuple<Rectangle, double> max_rectangle(double m_a, double b_b)  {
            //Initialize with list of maximum intervals.
            assert(!slab_heap.empty());
            auto& curr_splits = slab_heap[0].split_offsets;
            assert(!curr_splits.empty());

            std::vector<MaxIntervalAlt> initial_intervals;
            initial_intervals.reserve(curr_splits.size() - 1);

            for (size_t i = 0; i < curr_splits.size() - 1; ++i) {
                initial_intervals.emplace_back(curr_splits[i], curr_splits[i + 1]);
            }

            std::vector<std::tuple<Slab*, Slab*, Slab*, Slab*, pt2_t*, pt2_t*, std::vector<MaxIntervalAlt>>> slab_stack;

            assert(slab_heap[0].up != nullptr && slab_heap[0].down != nullptr);

            slab_stack.emplace_back(slab_heap[0].up->up, slab_heap[0].up->down,
                    slab_heap[0].down->up, slab_heap[0].down->down,
                    slab_heap[0].up->top_y, slab_heap[0].down->bottom_y,
                    initial_intervals
                    );


            auto update_weight = [](std::vector<MaxIntervalAlt>& max_intervals, std::vector<pt2_t*> merges, std::vector<double> weights, double scale) {
                size_t j = 0;
                for (size_t i = 0; i < merges.size(); ++i) {
                    while (j < max_intervals.size() && merges[i] != max_intervals[j].right()) {
                        j++;
                    }
                    max_intervals[j] += MaxIntervalAlt(merges[i], scale * weights[i]);
                }
            };

            auto merge_splits = [](std::vector<MaxIntervalAlt> const& max_intervals, std::vector<pt2_t*> const& splits) {
                size_t j = 0;
                std::vector<MaxIntervalAlt> new_intervals;
                for (size_t i = 0; i < splits.size(); ++i) {
                    while (j < max_intervals.size() - 1 && splits[i] != max_intervals[j].right()) {
                        new_intervals.emplace_back(max_intervals[j]);
                        j++;
                    }
                    new_intervals.emplace_back(max_intervals[j]);
                    new_intervals[j] += max_intervals[j + 1];
                    j++;
                }
                return new_intervals;
            };

            auto update_and_merge = [&](Slab* top, Slab* bottom, std::vector<MaxIntervalAlt>& m1) {
                auto m2 = m1;
                if (bottom != nullptr) { //if bottom is nullptr then we merge top or don't merge top.
                    update_weight(m1, bottom->m_merges, bottom->m_weights, m_a);
                    update_weight(m1, bottom->b_merges, bottom->b_weights, b_b);
                    m1 = merge_splits(m1, bottom->split_offsets);

                    m2 = merge_splits(m2, bottom->m_merges);
                    m2 = merge_splits(m2, bottom->b_merges);
                }
                if (top != nullptr) {
                    m1 = merge_splits(m1, top->m_merges);
                    m1 = merge_splits(m1, top->b_merges);
                    m2 = merge_splits(m2, top->split_offsets);
                }
                return m2;
            };
            auto get_top = [](Slab* p) {
                return p == nullptr ? nullptr : p->up;
            };
            auto get_bottom = [](Slab* p) {
                return p == nullptr ? nullptr : p->down;
            };

            auto get_mid = [](Slab* top, Slab* bottom) {
                //This assumes that both bottom and top cannot be null.
                return bottom == nullptr ? top->bottom_y : bottom->top_y;
            };

            double max_v = 0;
            Rectangle max_rect;
            while (!slab_stack.empty()) {
                auto [top_top, top_bottom, bottom_top, bottom_bottom, p_up, p_low, max_intervals] = slab_stack.back();
                slab_stack.pop_back();
                if (max_intervals.size() > 1) {
                    // m4 doesn't have any merges
                    // m2 is bottom_top merged.
                    // m3 is top_bottom merged.
                    // max_intervals is top_bottom and bottom_top merged.
                    assert(!(top_top == nullptr &&
                            top_bottom == nullptr &&
                            bottom_top == nullptr &&
                            bottom_bottom == nullptr));

                    if (!(top_top == nullptr && top_bottom == nullptr)) {
                        auto m2 = update_and_merge(top_top, top_bottom, max_intervals);
                        if (!(bottom_bottom == nullptr && bottom_top == nullptr)) {
                            auto m4 = update_and_merge(bottom_bottom, bottom_top, m2);
                            slab_stack.emplace_back(get_top(top_bottom), get_bottom(top_bottom),
                                    get_top(bottom_top), get_bottom(bottom_top),
                                    get_mid(top_top, top_bottom), get_mid(bottom_top, bottom_bottom), m4);
                        }
                        slab_stack.emplace_back(get_top(top_bottom), get_bottom(top_bottom),
                                                get_top(bottom_bottom), get_bottom(bottom_bottom),
                                                get_mid(top_top, top_bottom), p_low, m2);
                    }
                    if (!(bottom_bottom == nullptr && bottom_top == nullptr)) {
                        auto m3 = update_and_merge(bottom_bottom, bottom_top, max_intervals);
                        slab_stack.emplace_back(get_top(top_top), get_bottom(top_top),
                                                get_top(bottom_top), get_bottom(bottom_top),
                                                p_up, get_mid(bottom_top, bottom_bottom), m3);
                    }
                    slab_stack.emplace_back(get_top(top_top), get_bottom(top_top),
                                            get_top(bottom_bottom), get_bottom(bottom_bottom),
                                            p_up, p_low, max_intervals);

                }  else {
                    if (max_v < max_intervals[0].get_max().get_v()) {
                        double lx = max_intervals[0].get_max().get_l()->operator()(0);
                        double rx = max_intervals[0].get_max().get_r()->operator()(0);
                        double uy = (*p_up)(1);
                        double ly = (*p_low)(1);
                        max_rect = Rectangle(rx, uy, lx, ly);
                        max_v = max_intervals[0].get_max().get_v();
                    }
                }
            }
            return std::make_tuple(max_rect, max_v);

        }
    };


    std::tuple<Rectangle, double> max_rectangle(const wpoint_list_t& m_points, const wpoint_list_t& b_points, double eps, double a, double b) {
        SlabTreeAlt tree(m_points, b_points, eps);
        return tree.max_rectangle(a, b);
    }




}
