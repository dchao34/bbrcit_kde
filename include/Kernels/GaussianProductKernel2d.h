#ifndef BBRCITKDE_GAUSSIANPRODUCTKERNEL2D_H__
#define  BBRCITKDE_GAUSSIANPRODUCTKERNEL2D_H__

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
// + C(hx,hy) = 1 / (2 \pi * hx * hy)
//
// + U(p,hx,hy) = exp{-0.5 (px*px/(hx*hx) + py*py/(hy*hy)) }. 
//
// + U(p,q,hx,hy) := U(p-q,hx,hy)
//
template<typename T=double>
class CUDA_ALIGN16 GaussianProductKernel2d {

  public:

    using FloatType = T;

  public:
    CUDA_CALLABLE GaussianProductKernel2d();
    CUDA_CALLABLE GaussianProductKernel2d(T hx, T hy);

    CUDA_CALLABLE ~GaussianProductKernel2d() = default;
    CUDA_CALLABLE GaussianProductKernel2d(const GaussianProductKernel2d<T>&) = default;
    CUDA_CALLABLE GaussianProductKernel2d(GaussianProductKernel2d<T>&&) = default;
    CUDA_CALLABLE GaussianProductKernel2d& operator=(const GaussianProductKernel2d<T>&) = default;
    CUDA_CALLABLE GaussianProductKernel2d& operator=(GaussianProductKernel2d<T>&&) = default;

    // evaluate C(h,D)
    CUDA_CALLABLE T normalization() const;

    // evaluate U(x,y,h)
    template<typename PointT>
    CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&) const;

    // get/set the bandwidth
    CUDA_CALLABLE T hx() const;
    CUDA_CALLABLE T hy() const;
    CUDA_CALLABLE void set_hx(T);
    CUDA_CALLABLE void set_hy(T);
    CUDA_CALLABLE void set_bandwidths(T,T);

  private:
    T hx_;
    T hy_;

    // point_arg_eval: evaluates (px-qx)*(px-qx)/(hx*hx)+(py-qy)*(py-qy)/(hy*hy)
    // default behavior is provided through the function template, while 
    // specialized behavior are provided through overloads. 
    template<typename PointT>
      T point_arg_eval(const PointT&, const PointT&) const;

    CUDA_CALLABLE T point_arg_eval(const Point2d<T>&, const Point2d<T>&) const;
    
};

// Implementations
// ---------------

template<typename T>
GaussianProductKernel2d<T>::GaussianProductKernel2d() : hx_(1.0), hy_(1.0) {}

template<typename T>
GaussianProductKernel2d<T>::GaussianProductKernel2d(T hx, T hy) : hx_(hx), hy_(hy) {}

template<typename T>
inline void GaussianProductKernel2d<T>::set_hx(T hx) { hx_ = hx; }

template<typename T>
inline void GaussianProductKernel2d<T>::set_hy(T hy) { hy_ = hy; }

template<typename T>
inline void GaussianProductKernel2d<T>::set_bandwidths(T hx, T hy) {
  hx_ = hx; hy_ = hy;
}

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<typename T>
  template<typename PointT>
inline T GaussianProductKernel2d<T>::unnormalized_eval(
    const PointT &p, const PointT &q) const {
  return GaussianProduct2dTraits<T>::kernel(point_arg_eval(p, q));
}

template <typename T>
inline T GaussianProductKernel2d<T>::normalization() const {
  return GaussianProduct2dTraits<T>::normalization(hx_, hy_);
}

template<typename T>
  template<typename PointT>
inline T GaussianProductKernel2d<T>::point_arg_eval(
    const PointT &lhs, const PointT &rhs) const {
  return (lhs[0]-rhs[0])*(lhs[0]-rhs[0])/(hx_*hx_) +
         (lhs[1]-rhs[1])*(lhs[1]-rhs[1])/(hy_*hy_);
}

template<typename T>
inline T GaussianProductKernel2d<T>::point_arg_eval(
    const Point2d<T> &lhs, const Point2d<T> &rhs) const {
  return (lhs.x()-rhs.x())*(lhs.x()-rhs.x())/(hx_*hx_) +
         (lhs.y()-rhs.y())*(lhs.y()-rhs.y())/(hy_*hy_);
}

}

#endif
