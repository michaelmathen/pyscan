//
// Created by mmath on 7/7/17.
//
#include <random>

#include "Utilities.hpp"

#include "BloomFilter.hpp"

namespace pyscan {

  #define FNV_32_OFFSET 0x811c9dc5
  #define FNV_32_PRIME 16777619

  uint32_t fnv32a (uint32_t index) {
      uint32_t hash = FNV_32_OFFSET;
      hash = (hash ^ ((index >> 0) & 0xff)) * FNV_32_PRIME;
      hash = (hash ^ ((index >> 8) & 0xff)) * FNV_32_PRIME;
      hash = (hash ^ ((index >> 16) & 0xff)) * FNV_32_PRIME;
      hash = (hash ^ ((index >> 24) & 0xff)) * FNV_32_PRIME;
      return hash;
  }

  uint32_t fnv32a(uint32_t index, uint32_t seed) {
      uint32_t hash = seed;
      hash = (hash ^ ((index >> 0) & 0xff)) * FNV_32_PRIME;
      hash = (hash ^ ((index >> 8) & 0xff)) * FNV_32_PRIME;
      hash = (hash ^ ((index >> 16) & 0xff)) * FNV_32_PRIME;
      hash = (hash ^ ((index >> 24) & 0xff)) * FNV_32_PRIME;
      return hash;
  }

  BloomFilter::BloomFilter(int n, double p) :
          k(static_cast<int>(ceil(- log(p) / log(2)) + .5)),
          seeds(k, 0),
          size(static_cast<int>(ceil(-n * log(p)/ (log(2)*log(2))) +.5)),
          bit_arrays(size, false) {
      //std::random_device rd;

      for (int i = 0; i < k; i++) {
          seeds[i] = fnv32a(rand());
      }
  }

  BloomFilter::BloomFilter(int n, double p, int seed) : BloomFilter(n, p) {
      uint32_t seed_last = fnv32a(seed);
      for (int i = 0; i < k; i++) {
          seeds[i] = fnv32a(seed_last, seed);
          seed_last = seeds[i];
      }
  }

  BloomFilter::BloomFilter() : k(0), seeds(), size(0), bit_arrays() {}

  void BloomFilter::insert(uint32_t index) {
      for (int i = 0; i < k; i++) {
          uint32_t hash = fnv32a(index, seeds[i]) % size;
          bit_arrays[hash] = true;
      }
  }

  bool BloomFilter::mightBePresent(uint32_t index) const {
      if (bit_arrays.size() == 0) {
          return false;
      }
      for (int i = 0; i < k; i++) {
          uint32_t hash = fnv32a(index, seeds[i]) % size;
          if (!bit_arrays[hash]) { // any of the indices are 0 then we definitely haven't found it
              return false;
          }
      }
      return true; // Could be present unless we are very unlucky.
  }


  void BloomFilter::print(std::ostream& os) const {
      os << "BloomFilter " << bit_arrays;
  }
}
