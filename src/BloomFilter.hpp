//
// Created by mmath on 7/7/17.
//

#ifndef PYSCAN_BLOOMFILTER_HPP
#define PYSCAN_BLOOMFILTER_HPP
#include <cstdint>
#include <vector>

namespace pyscan {
  class BloomFilter {
      int k;
      std::vector<uint32_t> seeds;
      int size;
      std::vector<bool> bit_arrays;
  public:
      BloomFilter(int n, double p);
      BloomFilter(int n, double p, int seed);
      BloomFilter();

      void insert(uint32_t index);
      bool mightBePresent(uint32_t index) const;

      void print(std::ostream& os) const;
  };
}

#endif //PYSCAN_BLOOMFILTER_HPP
