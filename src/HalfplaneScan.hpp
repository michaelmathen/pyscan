//
// Created by mmath on 9/19/17.
//

#ifndef PYSCAN_HALFPLANE_HPP
#define PYSCAN_HALFPLANE_HPP

#include <functional>
#include <tuple>
#include "Point.hpp"

namespace pyscan {

    Pt3 lift_halfspace(Pt2 const& h, Pt3 const& p);
    Pt2 dropPoint(Pt3 const& fixed_point, Pt3 const& p1);

    std::tuple<Pt2, double> max_halfplane(
            point_list point_net,
            point_list red,
            weight_list red_w,
            point_list blue,
            weight_list blue_w,
            std::function<double(double, double)> const& f);


    std::tuple<Pt2, double> max_halfplane_simple(
            point_list point_net,
            point_list red,
            weight_list red_w,
            point_list blue,
            weight_list blue_w,
            std::function<double(double, double)> const& f);  

    std::tuple<Pt3, double> max_halfspace(
            point3_list point_net,
            point3_list red,
            weight_list red_w,
            point3_list blue,
            weight_list blue_w,
            std::function<double(double, double)> const& f);

    double under_line(Pt2& line, point_list& points, weight_list& weights);
    
    double evaluate_line(Pt2 const& line,
                        point_list const& red,
                        weight_list const& r_weight,
                        point_list const& blue,
                        weight_list const& b_weights,
                        std::function<double(double, double)> const& f);

  double evaluate_line_labeled(Pt2 const& line,
                        point_list const& red,
                        weight_list const& red_w,
                        label_list const& red_labels,
                        point_list const& blue,
                        weight_list const& blue_w,
                        label_list const& blue_labels,
                        std::function<double(double, double)> const& f);


   std::tuple<Point<>, double> max_halfplane_labeled(
            point_list point_net,
            point_list red,
            weight_list red_w,
            label_list red_labels,
            point_list blue,
            weight_list blue_w,
            label_list blue_labels,
            std::function<double(double, double)> const& f);

    std::tuple<Pt2, double> max_halfplane_simple_labeled(
            point_list point_net,
            point_list red,
            weight_list red_w,
            label_list red_labels,
            point_list blue,
            weight_list blue_w,
            label_list blue_labels,
            std::function<double(double, double)> const& f);

    std::tuple<Pt3, double> max_halfspace_labeled(
            point3_list point_net,
            point3_list red,
            weight_list red_w,
            label_list red_labels,
            point3_list blue,
            weight_list blue_w,
            label_list blue_labels,
            std::function<double(double, double)> const& f);
}

#endif //PYSCAN_HALFPLANE_HPP
