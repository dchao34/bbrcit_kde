#ifndef BBRCITKDE_GAUSSIANCONVKERNEL1D_H__
#define  BBRCITKDE_GAUSSIANCONVKERNEL1D_H__

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

// The one/two-point Gaussian convolution kernel in 1 dimension 
// with bandwidth hx is defined as follows:
//
// K(p;hx,hy) = C(hx) * U(p,hx)
// K(p,q;hx) := C(hx) * U(p,q,hx)
//
// where: 
//
// + C(hx) = 1 / (sqrt(4pi) / hx
//
// + U(p,hx,hy) = exp(-x*x/(4*hx*hx))
//
// + U(p,q,hx,hy) := U(p-q,hx,hy)
//
template<typename T=double>
class CUDA_ALIGN8 GaussianConvKernel1d {

  public:

    using FloatType = T;

  public:
    CUDA_CALLABLE GaussianConvKernel1d();
    CUDA_CALLABLE GaussianConvKernel1d(T hx);

    CUDA_CALLABLE ~GaussianConvKernel1d() = default;
    CUDA_CALLABLE GaussianConvKernel1d(const GaussianConvKernel1d<T>&) = default;
    CUDA_CALLABLE GaussianConvKernel1d(GaussianConvKernel1d<T>&&) = default;
    CUDA_CALLABLE GaussianConvKernel1d& operator=(const GaussianConvKernel1d<T>&) = default;
    CUDA_CALLABLE GaussianConvKernel1d& operator=(GaussianConvKernel1d<T>&&) = default;

    // evaluate C(hx,hy)
    CUDA_CALLABLE T normalization() const;

    // evaluate U(p,q,hx,hy)
    template<typename PointT>
      CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&) const;

    // evaluate U(p,q,hx*a,hy*a)
    template<typename PointT>
      CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&, T) const;

    // get/set the bandwidth
    CUDA_CALLABLE T hx() const;
    CUDA_CALLABLE void set_hx(T);

  private:
    T hx_;

    // point_arg_eval: evaluates (x-y)'(x-y)/(h*h*a*a).
    // default behavior is provided through the function template, while
    // specialized behavior are provided through overloads.
    template<typename PointT>
      T point_arg_eval(const PointT&, const PointT&, T a) const;
    CUDA_CALLABLE T point_arg_eval(const Point1d<T>&, const Point1d<T>&, T a) const;
    
};

// Implementations
// ---------------

template<typename T>
GaussianConvKernel1d<T>::GaussianConvKernel1d() : hx_(1.0) {}

template<typename T>
GaussianConvKernel1d<T>::GaussianConvKernel1d(T hx) : hx_(hx) {}

template<typename T>
inline void GaussianConvKernel1d<T>::set_hx(T hx) { hx_ = hx; }

template <typename T>
inline T GaussianConvKernel1d<T>::normalization() const {
  return GaussianConv1dTraits<T>::normalization(hx_);
}

template<typename T>
  template<typename PointT>
inline T GaussianConvKernel1d<T>::unnormalized_eval(
    const PointT &p, const PointT &q) const {
  return unnormalized_eval(p, q, ConstantTraits<T>::one());
}

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<typename T>
  template<typename PointT>
inline T GaussianConvKernel1d<T>::unnormalized_eval(
    const PointT &p, const PointT &q, T a) const {
  return GaussianConv1dTraits<T>::kernel(point_arg_eval(p,q,a));
}

template<typename T>
  template<typename PointT>
inline T GaussianConvKernel1d<T>::point_arg_eval(
    const PointT &lhs, const PointT &rhs, T a) const {
  return (lhs[0]-rhs[0])*(lhs[0]-rhs[0]) / (hx_*hx_*a*a);
}

template<typename T>
inline T GaussianConvKernel1d<T>::point_arg_eval(
    const Point1d<T> &lhs, const Point1d<T> &rhs, T a) const {
  return (lhs.x()-rhs.x())*(lhs.x()-rhs.x()) / (hx_*hx_*a*a);
}

}

#endif
