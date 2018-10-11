//
// Created by mmath on 5/28/17.
//
#include <functional>
#include <tuple>
#include <unordered_set>

#include "Range.hpp"
#include "RectangleScan.hpp"
#include "FunctionApprox.hpp"

namespace pyscan {


    /*
     * Takes a grid resolution,r, and a red set of points and a blue set of points. Computes a grid from this.z
     */
    Grid::Grid(size_t r_arg,
        point_list_t& red_points,
        weight_list_t& red_weight,
        point_list_t& blue_points,
        weight_list_t& blue_weight) :
        r(r_arg),
        red_counts(r_arg * r_arg, 0),
        blue_counts(r_arg * r_arg, 0),
        x_coords(),
        y_coords()
    {

        //Compute the grid from the n_begin and n_end and then fill in the values with the two sampled things.
        auto compX = [](Point<> const &p1, Point<> const &p2) {
            return p1(0) < p2(0);
        };
        auto compY = [](Point<> const &p1, Point<> const &p2) {
            return p1(1) < p2(1);
        };


        auto r_begin = red_points.begin();
        auto r_end = red_points.end();

        auto b_begin = blue_points.begin();
        auto b_end = blue_points.end();

        std::vector<pt2_t> x_pts(r, pt2_t());
        std::vector<pt2_t> y_pts(r, pt2_t());
        util::quantiles(r_begin, r_end, x_pts.begin(), x_pts.begin() + r / 2, red_weight.begin(), compX);
        util::quantiles(b_begin, b_end, x_pts.begin() + r / 2, x_pts.end() , blue_weight.begin(), compX);
        util::quantiles(r_begin, r_end, y_pts.begin(), y_pts.begin() + r / 2, red_weight.begin(), compY);
        util::quantiles(b_begin, b_end, x_pts.begin() + r / 2, x_pts.end() , blue_weight.begin(), compY);

        for (auto &el : x_pts) {
            x_coords.push_back(el(0));
        }
        for (auto &el : y_pts) {
            y_coords.push_back(el(1));
        }
        std::sort(x_coords.begin(), x_coords.end());
        std::sort(y_coords.begin(), y_coords.end());
        auto rw_it = red_weight.begin();
        for (auto point_it = r_begin; point_it != r_end; point_it++, rw_it++) {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), (*point_it)(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), (*point_it)(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                red_counts[iy * r + ix] += *rw_it;
            }
            total_red_weight += *rw_it;
        }
        auto bw_it = blue_weight.begin();
        for (auto point_it = b_begin; point_it != b_end; point_it++) {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), (*point_it)(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), (*point_it)(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                blue_counts[iy * r + ix] += *bw_it;
            }
            total_blue_weight += *bw_it;
        }
    }

    Grid::Grid(point_list_t& net, point_list_t& red, weight_list_t& red_w, point_list_t& blue, weight_list_t& blue_w) :
            r(net.size()),
            red_counts(r * r, 0),
            blue_counts(r * r, 0),
            x_coords(),
            y_coords() {

        for_each(net.begin(), net.end(), [&](Point<> const& pt) {
            x_coords.push_back((pt)(0));
        });
        for_each(net.begin(), net.end(), [&] (Point<> const& pt) {
            y_coords.push_back((pt)(1));
        });
        std::sort(x_coords.begin(), x_coords.end());
        std::sort(y_coords.begin(), y_coords.end());
        auto r_w_it = red_w.begin();
        for (auto& pt : red) {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), pt(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), pt(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                red_counts[iy * r + ix] += *r_w_it;
            }
            total_red_weight += *r_w_it;
            r_w_it++;
        }
        auto b_w_it = blue_w.begin();
        for (auto& pt : blue) {
            long ix = std::upper_bound(x_coords.begin(), x_coords.end(), pt(0)) - x_coords.begin() - 1;
            long iy = std::upper_bound(y_coords.begin(), y_coords.end(), pt(1)) - y_coords.begin() - 1;
            if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                blue_counts[iy * r + ix] += *b_w_it;
            }
            total_blue_weight += *b_w_it;
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
                         xCoord(sg.lowX()), yCoord(sg.lowY()),
                         sg.fValue());
    }

    /*
     * Simple 1/eps^4 algorithm described in the paper that just computes every subgrid.
     * This will work on a nonlinear function.
     */
    Subgrid maxSubgridNonLinear(Grid const &grid, std::function<double(double, double)> const& func) {
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
                        double maxf = func(red_l_count / t_red, blue_l_count / t_blue);
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


    Subgrid maxSubgridLinearSimple(Grid const &grid, double a, double b) {
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

    Subgrid maxSubgridLinearSimple(Grid const& grid, double eps,
        std::function<double(double, double)> const& f) {
        /*
         * This uses the
         */
        auto phi = [&] (Vec2 const& v) {return f(v[0], v[1]); };
        Subgrid curr_subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
        Subgrid max_subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
        double maxV = 0;
        auto linemaxF = [&] (Vec2 const& dir) {
            curr_subgrid = maxSubgridLinearSimple(grid, dir[0], dir[1]);
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


    Subgrid maxSubgridLinearTheory(Grid const& grid, double eps, 
        std::function<double(double, double)> const& f) {
        /*
         * This uses the
         */
        auto phi = [&] (Vec2 const& v) {return f(v[0], v[1]); };
        Subgrid curr_subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
        Subgrid max_subgrid(0, 0, 0, 0, -std::numeric_limits<double>::infinity());
        double maxV = 0;
        auto linemaxF = [&] (Vec2 const& dir) {
            curr_subgrid = maxSubgridLinearG(grid, grid.size() * r_ratio(grid.size()), dir[0], dir[1]);
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

    class LabelW {
       double weight;
       size_t label;

    public:
        size_t get_label() {
            return label;
        }

        double get_weight() {
            return weight;
        }


    };




    class LabeledGrid {

        std::vector<double> x_values;
        std::vector<double> y_values;
        std::vector<std::vector<LabelW>> labels;


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
                std::vector<double> partitions;
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



            auto part_m_x = accum_and_partition(m_points, m_total / r, 0);
            auto part_b_x = accum_and_partition(b_points, b_total / r, 0);
            x_values.reserve()
            std::merge(part_m_x.begin(), )
            std::sort(m_points.begin(), m_points.end(), compX);




            for (auto &el : x_pts) {
                x_coords.push_back(el(0));
            }
            for (auto &el : y_pts) {
                y_coords.push_back(el(1));
            }
            std::sort(x_coords.begin(), x_coords.end());
            std::sort(y_coords.begin(), y_coords.end());
            auto rw_it = red_weight.begin();
            for (auto point_it = r_begin; point_it != r_end; point_it++, rw_it++) {
                long ix = std::upper_bound(x_coords.begin(), x_coords.end(), (*point_it)(0)) - x_coords.begin() - 1;
                long iy = std::upper_bound(y_coords.begin(), y_coords.end(), (*point_it)(1)) - y_coords.begin() - 1;
                if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                    red_counts[iy * r + ix] += *rw_it;
                }
                total_red_weight += *rw_it;
            }
            auto bw_it = blue_weight.begin();
            for (auto point_it = b_begin; point_it != b_end; point_it++) {
                long ix = std::upper_bound(x_coords.begin(), x_coords.end(), (*point_it)(0)) - x_coords.begin() - 1;
                long iy = std::upper_bound(y_coords.begin(), y_coords.end(), (*point_it)(1)) - y_coords.begin() - 1;
                if (ix < r && iy < r && 0 <= ix && 0 <= iy) {
                    blue_counts[iy * r + ix] += *bw_it;
                }
                total_blue_weight += *bw_it;
            }

        }


    };


     Rectangle max_rect_labels(lpoint_list_t m_points, lpoint_list_t b_points, const discrepancy_func_t& func) {

         double m_Total = computeTotal(m_points);
         double b_Total = computeTotal(b_points);

         std::sort(m_points.begin(), m_points.end(), [](const pt2_t& p1, const pt2_t& p2){
            return p1(0) < p2(0);
         });

         std::sort(m_points.begin(), m_points.end(), [](const pt2_t& p1, const pt2_t& p2){
             return p1(0) < p2(0);
         });


         Rectangle maxRect(0, 0, 0, 0, 0.0);
         //Create a top and bottom line
         for (auto onb1 = net.begin(); onb1 != (net.end() - 2); onb1++) {
             for (auto onb2 = onb1 + 1; onb2 != (net.end() - 1); onb2++) {
                 auto nb1 = get<1>(*onb1) >= get<1>(*onb2) ? onb1 : onb2;
                 auto nb2 = get<1>(*onb1) < get<1>(*onb2) ? onb1 : onb2;

                 std::vector<LPoint<>> divisions;
                 std::vector<LPoint<>> m_slab;
                 std::vector<LPoint<>> b_slab;
                 auto between_f = [&](LPoint<> const& pt){
                     return get<1>(pt) <= get<1>(*nb1) && get<1>(*nb2) <= get<1>(pt);
                 };
                 auto pointXComp = [](LPoint<> const& p1, LPoint<> const& p2){
                     return get<0>(p1) < get<0>(p2);
                 };
                 std::remove_copy_if(net.begin(), net.end(), std::back_inserter(divisions), between_f);
                 //Filter the m_slab and b_slab points
                 std::remove_copy_if(m_points.begin(), m_points.end(), std::back_inserter(m_slab), between_f);
                 std::remove_copy_if(b_points.begin(), b_points.end(), std::back_inserter(b_slab), between_f);
                 std::sort(divisions.begin(), divisions.end(), pointXComp);

                 //Now position the points into each division
                 using part_t = std::vector<std::vector<LPoint<>>>;
                 part_t m_parts(divisions.size() + 1,  std::vector<LPoint<>>());
                 part_t b_parts(divisions.size() + 1, std::vector<LPoint<>>());

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

     Rectangle maxRectStatLabels(std::vector<LPoint<>> const& net,
                                  std::vector<LPoint<>> const& m_points,
                                  std::vector<LPoint<>> const& b_points,
                                  double rho) {
         return maxRectLabels(net, m_points, b_points, [&](double mr, double br){
             return kulldorff(mr, br, rho);
         });
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


    /*
    auto operator<<(std::ostream& os, const I_Type & t) -> std::ostream& {
        if (t == I_Type::VALUE) {
            os << "VALUE";
        } else {
            os << "MAX";
        }
        return os;
    }
     */
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
            stack_nodes.push_back(std::make_tuple(lower, upper, root));
            while (stack_nodes.size() > 0) {
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
                stack_nodes.push_back(std::make_tuple(std::get<0>(node), mid, left_child));
                stack_nodes.push_back(std::make_tuple(mid + 1, std::get<1>(node), right_child));
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

        while (!sub_problem_stack.empty()) {
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
            while (!slabs.empty()) {
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

    Subgrid maxSubgridLinearG(Grid const &grid, long r_prime, double a, double b) {
        SlabTree slabT(grid, (size_t)r_prime);
#ifdef DEBUG
        std::cout<<"finished the slab tree" << std::endl;
#endif
        //std::cout << slabT << std::endl;
        MaximumIntervals maxInt(grid.size());
        return maxSubgridLinear(slabT.root, maxInt, a, b);
    }


}
