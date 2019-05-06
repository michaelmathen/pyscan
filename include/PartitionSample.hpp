//
// Created by mmath on 12/9/18.
//

#ifndef PYSCAN_PARTITIONSAMPLE_HPP
#define PYSCAN_PARTITIONSAMPLE_HPP
#include "Point.hpp"

namespace pyscan {
    //Gets a line that is within .1 of the max line
    static const size_t HAM_NET_SIZE = 10;

    wpoint_list_t ham_tree_sample(const wpoint_list_t &pts, size_t s_size);

   // lpoint_list_t ham_tree_sample(const lpoint_list_t &pts, size_t s_size);
}
#endif //PYSCAN_PARTITIONSAMPLE_HPP
