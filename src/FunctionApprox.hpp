//
// Created by mmath on 7/7/17.
//

#ifndef PYSCAN_FUNCTIONAPPROX_HPP
#define PYSCAN_FUNCTIONAPPROX_HPP
#include <functional>
#include "Vecky.hpp"

namespace pyscan {

  double approximateHull(double eps,
                        VecD const& cc, VecD const& cl,
                        std::function<double(VecD)> phi, //function to maximize
                        std::function<VecD(VecD)> lineMaxF);

  double approximateHull(double eps,
            std::function<double(VecD)> phi, //function to maximize
            std::function<VecD(VecD)> lineMaxF);
}
#endif //PYSCAN_FUNCTIONAPPROX_HPP
