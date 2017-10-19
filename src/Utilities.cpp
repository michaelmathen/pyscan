#include <cstdint>

#include "Utilities.hpp"

namespace pyscan {
    double invsqrt( double number ) {
        double y = number;
        double x2 = y * 0.5;
        std::int64_t i = *(std::int64_t *) &y;
        // The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
        i = 0x5fe6eb50c7b537a9 - (i >> 1);
        y = *(double *) &i;
        y = y * (1.5 - (x2 * y * y));   // 1st iteration
        y  = y * ( 1.5 - ( x2 * y * y ) );   // 2nd iteration, this can be removed
        return y;
    }



}
