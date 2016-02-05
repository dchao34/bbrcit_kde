#ifndef BBRCITKDE_FLOATUTILS_H__
#define BBRCITKDE_FLOATUTILS_H__

#include <cmath>
#include <limits>
#include <type_traits>

namespace bbrcit {

// Idea taken from:
// http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template <typename T> 
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type 
cppref_almost_equal(T lhs, T rhs, int ulp) {
  // case 1: relative epsilon comparison for normalized numbers. 
  // case 2: absolute error comparison for subnormal numbers.
  return std::abs(lhs-rhs) < std::numeric_limits<T>::epsilon() * ulp * std::abs(lhs+rhs)
         || std::abs(lhs-rhs) < std::numeric_limits<T>::min();
}

// compares whether two floating numbers are ``nearly equal''. in particular,
// nearly_equal() returns true iff one of the following hold:
// 1. their absolute error is at most the minimum normalized number.
// 2. their relative error is at most the machine epsilon. 
//
// Idea taken from:
// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
//
// A similar idea: 
// http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template <typename T> 
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type 
almost_equal(T lhs, T rhs) {

  T diff = std::abs(lhs-rhs);

  // absolute error
  if (diff <= std::numeric_limits<T>::min()) { return true; }

  // relative error
  return diff <= std::numeric_limits<T>::epsilon() * std::max(std::abs(lhs), std::abs(rhs));
}

// compares whether two floating numbers are approximately equal up to
// some relative error and absolute error. 
// approximately_equal() returns true iff one of the following hold:
// 1. their absolute error is at most abs_err.
// 2. their relative error is at most rel_err. 
template <typename T> 
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type 
approximately_equal(T lhs, T rhs, 
                   T rel_err, T abs_err = std::numeric_limits<T>::min()) {

  T diff = std::abs(lhs-rhs);

  // absolute error
  if (diff <= abs_err) { return true; }

  // relative error
  return diff <= rel_err * std::max(std::abs(lhs), std::abs(rhs));
}

template<typename PointT> 
bool ExactEqual(const PointT &lhs, const PointT &rhs) {
  int i = 0; while (i < lhs.dim() && lhs[i] == rhs[i]) { ++i; }
  return i == lhs.dim();
}

// lexicographic comparison of PointT objects. note that the 
// equality comparison is == even for floats; this is intentional. 
template<typename PointT> 
bool ExactLexicoLess(const PointT &lhs, const PointT &rhs) {
  int i = 0; while (i < lhs.dim() && lhs[i] == rhs[i]) { ++i; }
  return i != lhs.dim() && lhs[i] < rhs[i];
}

}

#endif
