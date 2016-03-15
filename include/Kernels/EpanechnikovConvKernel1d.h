#ifndef BBRCITKDE_EPANECHNIKOVCONVKERNEL1D_H__
#define  BBRCITKDE_EPANECHNIKOVCONVKERNEL1D_H__

#ifdef __CUDACC__
#define CUDA_ALIGN8 __align__(8)
#define CUDA_ALIGN16 __align__(16)
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_ALIGN8
#define CUDA_ALIGN16
#define CUDA_CALLABLE
#endif

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>

#include <Kernels/KernelTraits.h>
#include <KdeTraits.h>
#include <Point1d.h>

namespace bbrcit {

// The one/two-point Epanechnikov convolution kernel in 1 dimension 
// with bandwidth hx is defined as follows:
//
// K(p;hx,hy) = C(hx) * U(p,hx)
// K(p,q;hx) := C(hx) * U(p,q,hx)
//
// where: 
//
// + C(hx) = (3 / 160) / hx
//
// + U(p,hx,hy) = 
//       (2-|px/hx|)^3 * (|px/hx|^2+6*|px/hx|+4) if 0 <= |px/hx| < 2
//       0                                       otherwise
//
// + U(p,q,hx,hy) := U(p-q,hx,hy)
//
template<typename T=double>
class CUDA_ALIGN8 EpanechnikovConvKernel1d {

  public:

    using FloatType = T;

  public:
    CUDA_CALLABLE EpanechnikovConvKernel1d();
    CUDA_CALLABLE EpanechnikovConvKernel1d(T hx);

    CUDA_CALLABLE ~EpanechnikovConvKernel1d() = default;
    CUDA_CALLABLE EpanechnikovConvKernel1d(const EpanechnikovConvKernel1d<T>&) = default;
    CUDA_CALLABLE EpanechnikovConvKernel1d(EpanechnikovConvKernel1d<T>&&) = default;
    CUDA_CALLABLE EpanechnikovConvKernel1d& operator=(const EpanechnikovConvKernel1d<T>&) = default;
    CUDA_CALLABLE EpanechnikovConvKernel1d& operator=(EpanechnikovConvKernel1d<T>&&) = default;

    // evaluate C(hx,hy)
    CUDA_CALLABLE T normalization() const;

    // evaluate U(p,q,hx,hy)
    template<typename PointT>
      CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&) const;

    // evaluate U(p,q,hx*a,hy*a)
    template<typename PointT>
      CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&, T) const;
    CUDA_CALLABLE T unnormalized_eval(const Point1d<T>&, const Point1d<T>&, T) const;

    // get/set the bandwidth
    CUDA_CALLABLE T hx() const;
    CUDA_CALLABLE void set_hx(T);

  private:
    T hx_;
    
};

// Implementations
// ---------------

template<typename T>
EpanechnikovConvKernel1d<T>::EpanechnikovConvKernel1d() : hx_(1.0) {}

template<typename T>
EpanechnikovConvKernel1d<T>::EpanechnikovConvKernel1d(T hx) : hx_(hx) {}

template<typename T>
inline void EpanechnikovConvKernel1d<T>::set_hx(T hx) { hx_ = hx; }

template <typename T>
inline T EpanechnikovConvKernel1d<T>::normalization() const {
  return EpanechnikovConv1dTraits<T>::normalization(hx_);
}

template<typename T>
  template<typename PointT>
inline T EpanechnikovConvKernel1d<T>::unnormalized_eval(
    const PointT &p, const PointT &q) const {
  return unnormalized_eval(p, q, ConstantTraits<T>::one());
}

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<typename T>
  template<typename PointT>
inline T EpanechnikovConvKernel1d<T>::unnormalized_eval(
    const PointT &lhs, const PointT &rhs, T a) const {
  return EpanechnikovConv1dTraits<T>::kernel((lhs[0]-rhs[0])/(hx_*a));
}

template<typename T>
inline T EpanechnikovConvKernel1d<T>::unnormalized_eval(
    const Point1d<T> &lhs, const Point1d<T> &rhs, T a) const {
  return EpanechnikovConv1dTraits<T>::kernel((lhs.x()-rhs.x())/(hx_*a));
}

}

#endif
