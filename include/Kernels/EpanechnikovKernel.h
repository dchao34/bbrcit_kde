#ifndef BBRCITKDE_EPANECHNIKOVKERNEL_H__
#define BBRCITKDE_EPANECHNIKOVKERNEL_H__

#define _USE_MATH_DEFINES

#ifdef __CUDACC__
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_CALLABLE
#endif

#include <cmath>
#include <iostream>

#include <Kernels/KernelTraits.h>
#include <KdeTraits.h>
#include <Point2d.h>

namespace bbrcit {

// The one/two-point EpanechnikovKernel kernel in D dimensions of 
// bandwidth h is defined as follows:
//
// K(x;h,D) = C(h,D) * U(x,h)
// K(x,y;h,D) := C(h,D) * U(x,y,h)
//
// where: 
//
// + C(h,D) = 0.5 * (D+2) / (unit_volume(D) * h^D)
//   unit_volume(D) = \pi^{D/2} * \Gamma(1+D/2), where
//
// + U(x,h) = 
//       1 - x'x/(h*h), if x'x < 1
//       0            , otherwise
//
// + U(x,y,h) := U(x-y,h)
//
template<int D, typename T=double>
class EpanechnikovKernel {

  public:

    using FloatType = T;
    static constexpr int dim() { return D; }

  public:

    CUDA_CALLABLE EpanechnikovKernel();
    CUDA_CALLABLE EpanechnikovKernel(T bandwidth);

    CUDA_CALLABLE ~EpanechnikovKernel() = default;
    CUDA_CALLABLE EpanechnikovKernel(const EpanechnikovKernel&) = default;
    CUDA_CALLABLE EpanechnikovKernel(EpanechnikovKernel&&) = default;
    CUDA_CALLABLE EpanechnikovKernel& operator=(const EpanechnikovKernel&) = default;
    CUDA_CALLABLE EpanechnikovKernel& operator=(EpanechnikovKernel&&) = default;

    // evaluate C(h,D)
    CUDA_CALLABLE T normalization() const;

    // evaluate U(x,y,h)
    template<typename PointT>
    CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&) const;

    // get/set the bandwidth
    CUDA_CALLABLE T bandwidth() const;
    CUDA_CALLABLE void set_bandwidth(T);

  private: 
    T bandwidth_;

    // point_arg_eval: evaluates (x-y)'(x-y)/h*h. 
    // default behavior is provided through the function template, while 
    // specialized behavior are provided through overloads. 
    template<typename PointT>
      T point_arg_eval(const PointT&, const PointT&) const;

    CUDA_CALLABLE T point_arg_eval(const Point2d<T>&, const Point2d<T>&) const;

    
};

// Implementations
// ---------------

template<int D, typename T>
EpanechnikovKernel<D,T>::EpanechnikovKernel() : bandwidth_(1.0) {}

template<int D, typename T>
EpanechnikovKernel<D,T>::EpanechnikovKernel(T bw) : bandwidth_(bw) {}

template<int D, typename T>
inline T EpanechnikovKernel<D,T>::bandwidth() const { return bandwidth_; }

template<int D, typename T>
inline void EpanechnikovKernel<D,T>::set_bandwidth(T bw) { bandwidth_ = bw; }

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<int D, typename T>
  template<typename PointT>
inline T EpanechnikovKernel<D,T>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return EpanechnikovTraits<D,T>::kernel(point_arg_eval(p, q));
}

template <int D, typename T>
inline T EpanechnikovKernel<D,T>::normalization() const {
  return EpanechnikovTraits<D,T>::normalization(bandwidth_);
}

template<int D, typename T>
  template<typename PointT>
T EpanechnikovKernel<D,T>::point_arg_eval(const PointT &lhs, const PointT &rhs) const {

  T result = ConstantTraits<T>::zero(); 

  T diff; 
  for (int i = 0; i < D; ++i) {
    diff = lhs[i] - rhs[i];
    result += diff * diff;
  }

  return result / (bandwidth_ * bandwidth_);
}

template<int D, typename T>
inline T EpanechnikovKernel<D,T>::point_arg_eval(
    const Point2d<T> &lhs, 
    const Point2d<T> &rhs) const {
  return ((lhs.x()-rhs.x())*(lhs.x()-rhs.x()) +
          (lhs.y()-rhs.y())*(lhs.y()-rhs.y())) / (bandwidth_*bandwidth_);
}


}

#endif
