//
// Created by mmath on 6/28/17.
//

#ifndef PYSCAN_UTILITIES_HPP
#define PYSCAN_UTILITIES_HPP
#include <vector>
#include <sstream>
#include <ostream>

namespace pyscan {


  double invsqrt(double number);

  
  template<class T>
  auto operator<<(std::ostream& os, const T& t) -> decltype(t.print(os), os) {
      t.print(os);
      return os;
  }


  template <typename T>
  std::ostream& operator<< (std::ostream& out, std::vector<T> const& els) {
      out << "[";
      for (auto i : els) {
          out << i << ", ";
      }
      out << "]";
      return out;
  }


}
#endif //PYSCAN_UTILITIES_HPP
