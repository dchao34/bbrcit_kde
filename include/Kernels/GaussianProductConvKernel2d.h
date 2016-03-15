#ifndef BBRCITKDE_GAUSSIANPRODUCTCONVKERNEL2D_H__
#define  BBRCITKDE_GAUSSIANPRODUCTCONVKERNEL2D_H__

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

// The one/two-point Gaussian product kernel in 2 dimensions 
// with bandwidth hx and hy is defined as follows:
//
// K(p;hx,hy) = C(hx,hy) * U(p,hx,hy)
// K(p,q;hx,hy) := C(hx,hy) * U(p,q,hx,hy)
//
// where: 
//
// + C(hx,hy) = 1 / (4 \pi * hx * hy)
//
// + U(p,hx,hy) = exp{-0.25 (px*px/(hx*hx) + py*py/(hy*hy)) }. 
//
// + U(p,q,hx,hy) := U(p-q,hx,hy)
//
template<typename T=double>
class CUDA_ALIGN16 GaussianProductConvKernel2d {

  public:

    using FloatType = T;

  public:
    CUDA_CALLABLE GaussianProductConvKernel2d();
    CUDA_CALLABLE GaussianProductConvKernel2d(T hx, T hy);

    CUDA_CALLABLE ~GaussianProductConvKernel2d() = default;
    CUDA_CALLABLE GaussianProductConvKernel2d(const GaussianProductConvKernel2d<T>&) = default;
    CUDA_CALLABLE GaussianProductConvKernel2d(GaussianProductConvKernel2d<T>&&) = default;
    CUDA_CALLABLE GaussianProductConvKernel2d& operator=(const GaussianProductConvKernel2d<T>&) = default;
    CUDA_CALLABLE GaussianProductConvKernel2d& operator=(GaussianProductConvKernel2d<T>&&) = default;

    // evaluate C(h,D)
    CUDA_CALLABLE T normalization() const;

    // evaluate U(x,y,h)
    template<typename PointT>
    CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&) const;

    // evaluate U(x,y,h*a)
    template<typename PointT>
    CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&, T) const;

    // get/set the bandwidth
    CUDA_CALLABLE T hx() const;
    CUDA_CALLABLE T hy() const;
    CUDA_CALLABLE void set_hx(T);
    CUDA_CALLABLE void set_hy(T);
    CUDA_CALLABLE void set_bandwidths(T,T);

  private:
    T hx_;
    T hy_;

    // point_arg_eval: evaluates (px-qx)*(px-qx)/(hx*hx*a*a)+(py-qy)*(py-qy)/(hy*hy*a*a)
    // default behavior is provided through the function template, while 
    // specialized behavior are provided through overloads. 
    template<typename PointT>
      T point_arg_eval(const PointT&, const PointT&, T) const;

    CUDA_CALLABLE T point_arg_eval(const Point2d<T>&, const Point2d<T>&, T) const;
    
};

// Implementations
// ---------------

template<typename T>
GaussianProductConvKernel2d<T>::GaussianProductConvKernel2d() : hx_(1.0), hy_(1.0) {}

template<typename T>
GaussianProductConvKernel2d<T>::GaussianProductConvKernel2d(T hx, T hy) : hx_(hx), hy_(hy) {}

template<typename T>
inline void GaussianProductConvKernel2d<T>::set_hx(T hx) { hx_ = hx; }

template<typename T>
inline void GaussianProductConvKernel2d<T>::set_hy(T hy) { hy_ = hy; }

template<typename T>
inline void GaussianProductConvKernel2d<T>::set_bandwidths(T hx, T hy) {
  hx_ = hx; hy_ = hy;
}

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<typename T>
  template<typename PointT>
inline T GaussianProductConvKernel2d<T>::unnormalized_eval(
    const PointT &p, const PointT &q) const {
  return unnormalized_eval(p, q, ConstantTraits<T>::one());
}

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<typename T>
  template<typename PointT>
inline T GaussianProductConvKernel2d<T>::unnormalized_eval(
    const PointT &p, const PointT &q, T a) const {
  return GaussianProductConv2dTraits<T>::kernel(point_arg_eval(p, q, a));
}

template <typename T>
inline T GaussianProductConvKernel2d<T>::normalization() const {
  return GaussianProductConv2dTraits<T>::normalization(hx_, hy_);
}

template<typename T>
  template<typename PointT>
inline T GaussianProductConvKernel2d<T>::point_arg_eval(
    const PointT &lhs, const PointT &rhs, T a) const {
  return (lhs[0]-rhs[0])*(lhs[0]-rhs[0])/(hx_*hx_*a*a) +
         (lhs[1]-rhs[1])*(lhs[1]-rhs[1])/(hy_*hy_*a*a);
}

template<typename T>
inline T GaussianProductConvKernel2d<T>::point_arg_eval(
    const Point2d<T> &lhs, const Point2d<T> &rhs, T a) const {
  return (lhs.x()-rhs.x())*(lhs.x()-rhs.x())/(hx_*hx_*a*a) +
         (lhs.y()-rhs.y())*(lhs.y()-rhs.y())/(hy_*hy_*a*a);
}

}

#endif
