#ifndef BBRCITKDE_GAUSSIANKERNEL_H__
#define  BBRCITKDE_GAUSSIANKERNEL_H__

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

// The one/two-point Gaussian kernel in D dimensions of 
// bandwidth h is defined as follows:
//
// K(x;h,D) = C(h,D) * U(x,h)
// K(x,y;h,D) := C(h,D) * U(x,y,h)
//
// where: 
//
// + C(h,D) = 1 / (unit_volume(D) * h^D)
//   unit_volume(D) = (2 \pi)^{-D/2}
//
// + U(x,h) = exp{-0.5 x'x/ (h*h) }. 
//   x'x := dot product of the D-dim vector x.
//
// + U(x,y,h) := U(x-y,h)
//
template<int D, typename T=double>
class GaussianKernel {

  public:

    using FloatType = T;
    static constexpr int dim() { return D; }

  public:
    CUDA_CALLABLE GaussianKernel();
    CUDA_CALLABLE GaussianKernel(T bandwidth);

    CUDA_CALLABLE ~GaussianKernel() = default;
    CUDA_CALLABLE GaussianKernel(const GaussianKernel&) = default;
    CUDA_CALLABLE GaussianKernel(GaussianKernel&&) = default;
    CUDA_CALLABLE GaussianKernel& operator=(const GaussianKernel&) = default;
    CUDA_CALLABLE GaussianKernel& operator=(GaussianKernel&&) = default;

    // evaluate C(h,D)
    CUDA_CALLABLE T normalization() const;

    // evaluate U(x,y,h)
    template<typename PointT>
    CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&) const;

    // evaluate U(x,y,h*a)
    template<typename PointT>
    CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&, T a) const;

    // get/set the bandwidth
    CUDA_CALLABLE T bandwidth() const;
    CUDA_CALLABLE void set_bandwidth(T);

  private:
    T bandwidth_;

    // point_arg_eval: evaluates (x-y)'(x-y)/(h*h*a*a). 
    // default behavior is provided through the function template, while 
    // specialized behavior are provided through overloads. 
    template<typename PointT>
      T point_arg_eval(const PointT&, const PointT&, T a) const;

    CUDA_CALLABLE T point_arg_eval(const Point2d<T>&, const Point2d<T>&, T a) const;
    
};

// Implementations
// ---------------

template<int D, typename T>
GaussianKernel<D, T>::GaussianKernel() : bandwidth_(1.0) {}

template<int D, typename T>
GaussianKernel<D,T>::GaussianKernel(T bw) : bandwidth_(bw) {}

template<int D, typename T>
inline T GaussianKernel<D,T>::bandwidth() const { return bandwidth_; }

template<int D, typename T>
inline void GaussianKernel<D,T>::set_bandwidth(T bw) { bandwidth_ = bw; }

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<int D, typename T>
  template<typename PointT>
inline T GaussianKernel<D,T>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return unnormalized_eval(p, q, ConstantTraits<T>::one());
}

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<int D, typename T>
  template<typename PointT>
inline T GaussianKernel<D,T>::unnormalized_eval(const PointT &p, const PointT &q, T a) const {
  return GaussianTraits<D,T>::kernel(point_arg_eval(p, q, a));
}

template <int D, typename T>
inline T GaussianKernel<D,T>::normalization() const {
  return GaussianTraits<D,T>::normalization(bandwidth_);
}

template<int D, typename T>
  template<typename PointT>
T GaussianKernel<D,T>::point_arg_eval(const PointT &lhs, const PointT &rhs, T a) const {

  T result = ConstantTraits<T>::zero(); 

  T diff; 
  for (int i = 0; i < D; ++i) {
    diff = lhs[i] - rhs[i];
    result += diff * diff;
  }

  return result / (bandwidth_ * bandwidth_ * a * a);
}

template<int D, typename T>
inline T GaussianKernel<D,T>::point_arg_eval(
    const Point2d<T> &lhs, const Point2d<T> &rhs, T a) const {
  return ((lhs.x()-rhs.x())*(lhs.x()-rhs.x()) +
          (lhs.y()-rhs.y())*(lhs.y()-rhs.y())) / (bandwidth_*bandwidth_*a*a);
}


}

#endif
