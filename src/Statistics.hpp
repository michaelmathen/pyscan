//
// Created by mmath on 7/7/17.
//

#ifndef PYSCAN_STATISTICS_HPP
#define PYSCAN_STATISTICS_HPP
#include <cmath>
#include <limits>

namespace pyscan {

    inline double kulldorff(double mr, double br, double rho) {
        if (mr < rho || br < rho || br > 1 - rho || mr > 1 - rho) {
            return 0;
        }
        if (abs(1 - abs(mr / br)) <= std::numeric_limits<double>::epsilon()) {
            return 0;
        }
        if (mr <= 0 || mr >= 1)
            return std::numeric_limits<double>::infinity();
        if (br <= 0 || br >= 1)
            return std::numeric_limits<double>::infinity();
        return mr * log(mr / br) + (1 - mr) * log((1 - mr) / (1 - br));
    }


    inline double gamma(double mr, double br, double rho) {
        if (br < rho || br > 1 - rho)
            return 0;
        if (br <= 0 || br >= 1)
            return std::numeric_limits<double>::infinity();
        return (mr - br) * (mr - br) / (br * (1 - br));
    }

};
#endif //PYSCAN_STATISTICS_HPP
