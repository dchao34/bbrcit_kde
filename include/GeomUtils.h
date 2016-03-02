#ifndef GEOM_UTILS_H__
#define GEOM_UTILS_H__

#include <KdeTraits.h>

namespace bbrcit {

// straight forward implementation of a dot product between two points. 
template<typename PointT, 
         typename FloatT = typename PointT::FloatType>
FloatT DotProduct(const PointT &lhs, const PointT &rhs) {

  FloatT result = ConstantTraits<FloatT>::zero();
  for (int i = 0; i < lhs.dim(); ++i) {
    result += lhs[i] * rhs[i];
  }

  return result;
}

}

#endif
