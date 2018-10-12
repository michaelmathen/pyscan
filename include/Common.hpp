#ifndef __COMMON_H__
#define __COMMON_H__

#include <tuple>
#include <functional>
#include <algorithm>
#include "Point.hpp"

namespace spscan{

using discrepancy_func_t = std::function<double(double, double)>;

}

#endif
