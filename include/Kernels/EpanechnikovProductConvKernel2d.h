#ifndef BBRCITKDE_EPANECHNIKOVPRODUCTCONVKERNEL2D_H__
#define  BBRCITKDE_EPANECHNIKOVPRODUCTCONVKERNEL2D_H__

#ifdef __CUDACC__
#define CUDA_ALIGN16 __align__(16)
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_ALIGN16
#define CUDA_CALLABLE
#endif

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>

#include <Kernels/KernelTraits.h>
#include <KdeTraits.h>
#include <Point2d.h>

namespace bbrcit {

// The one/two-point Epanechnikov product convolution kernel in 2 dimensions 
// with bandwidth hx and hy is defined as follows:
//
// K(p;hx,hy) = C(hx,hy) * U(p,hx,hy)
// K(p,q;hx,hy) := C(hx,hy) * U(p,q,hx,hy)
//
// where: 
//
// + C(hx,hy) = (3*3 / (160*160)) / (hx*hy)
//
// + U(p,hx,hy) = 
//       (2-px/hx)^3*((px/hx)^2+6*px/hx)+4)*
//         (2-py/hy)^3*((py/hy)^2+6*py/hy)+4), if 0<=|px/hx|<2 and 0<=|py/hy|<2 
//       0                                   , otherwise
//
// + U(p,q,hx,hy) := U(p-q,hx,hy)
//
template<typename T=double>
class CUDA_ALIGN16 EpanechnikovProductConvKernel2d {

  public:

    using FloatType = T;

  public:
    CUDA_CALLABLE EpanechnikovProductConvKernel2d();
    CUDA_CALLABLE EpanechnikovProductConvKernel2d(T hx, T hy);

    CUDA_CALLABLE ~EpanechnikovProductConvKernel2d() = default;
    CUDA_CALLABLE EpanechnikovProductConvKernel2d(const EpanechnikovProductConvKernel2d<T>&) = default;
    CUDA_CALLABLE EpanechnikovProductConvKernel2d(EpanechnikovProductConvKernel2d<T>&&) = default;
    CUDA_CALLABLE EpanechnikovProductConvKernel2d& operator=(const EpanechnikovProductConvKernel2d<T>&) = default;
    CUDA_CALLABLE EpanechnikovProductConvKernel2d& operator=(EpanechnikovProductConvKernel2d<T>&&) = default;

    // evaluate C(hx,hy)
    CUDA_CALLABLE T normalization() const;

    // evaluate U(p,q,hx,hy)
    template<typename PointT>
      CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&) const;

    // evaluate U(p,q,hx*a,hy*a)
    template<typename PointT>
      CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&, T) const;
    CUDA_CALLABLE T unnormalized_eval(const Point2d<T>&, const Point2d<T>&, T) const;

    // get/set the bandwidth
    CUDA_CALLABLE T hx() const;
    CUDA_CALLABLE T hy() const;
    CUDA_CALLABLE void set_hx(T);
    CUDA_CALLABLE void set_hy(T);
    CUDA_CALLABLE void set_bandwidths(T,T);

  private:
    T hx_;
    T hy_;
    
};

// Implementations
// ---------------

template<typename T>
EpanechnikovProductConvKernel2d<T>::EpanechnikovProductConvKernel2d() : hx_(1.0), hy_(1.0) {}

template<typename T>
EpanechnikovProductConvKernel2d<T>::EpanechnikovProductConvKernel2d(T hx, T hy) : hx_(hx), hy_(hy) {}

template<typename T>
inline void EpanechnikovProductConvKernel2d<T>::set_hx(T hx) { hx_ = hx; }

template<typename T>
inline void EpanechnikovProductConvKernel2d<T>::set_hy(T hy) { hy_ = hy; }

template<typename T>
inline void EpanechnikovProductConvKernel2d<T>::set_bandwidths(T hx, T hy) {
  hx_ = hx; hy_ = hy;
}

template <typename T>
inline T EpanechnikovProductConvKernel2d<T>::normalization() const {
  return EpanechnikovProductConv2dTraits<T>::normalization(hx_, hy_);
}

template<typename T>
  template<typename PointT>
inline T EpanechnikovProductConvKernel2d<T>::unnormalized_eval(
    const PointT &p, const PointT &q) const {
  return unnormalized_eval(p, q, ConstantTraits<T>::one());
}

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<typename T>
  template<typename PointT>
inline T EpanechnikovProductConvKernel2d<T>::unnormalized_eval(
    const PointT &lhs, const PointT &rhs, T a) const {
  return EpanechnikovProductConv2dTraits<T>::kernel((lhs[0]-rhs[0])/(hx_*a)) *
         EpanechnikovProductConv2dTraits<T>::kernel((lhs[1]-rhs[1])/(hy_*a));
}

template<typename T>
inline T EpanechnikovProductConvKernel2d<T>::unnormalized_eval(
    const Point2d<T> &lhs, const Point2d<T> &rhs, T a) const {
  return EpanechnikovProductConv2dTraits<T>::kernel((lhs.x()-rhs.x())/(hx_*a)) *
         EpanechnikovProductConv2dTraits<T>::kernel((lhs.y()-rhs.y())/(hy_*a));
}

}

#endif
