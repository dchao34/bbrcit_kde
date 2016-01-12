#include <cmath>
#include <limits>
#include <type_traits>

// Idea taken from:
// http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template <typename T> 
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type 
almost_equal(T lhs, T rhs, int ulp) {
  // case 1: relative epsilon comparison for normalized numbers. 
  // case 2: absolute error comparison for subnormal numbers.
  return std::abs(lhs-rhs) < std::numeric_limits<T>::epsilon() * ulp * std::abs(lhs+rhs)
         || std::abs(lhs-rhs) < std::numeric_limits<T>::min();
}
