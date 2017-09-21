//
// Created by mmath on 7/7/17.
//

#ifndef PYSCAN_FUNCTIONAPPROX_HPP
#define PYSCAN_FUNCTIONAPPROX_HPP
#include <functional>

double approximateHull(double mi, double bi,
                       double mj, double bj,
                       double alpha, double theta, double eps,
                       std::function<double(double, double)> phi, //function to maximize
                       std::function<std::tuple<double, double>(double)> lineMaxF);

#endif //PYSCAN_FUNCTIONAPPROX_HPP
