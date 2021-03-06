#ifndef BBRCITKDE_EPANECHNIKOVPRODUCTKERNEL2D_H__
#define  BBRCITKDE_EPANECHNIKOVPRODUCTKERNEL2D_H__

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
#include <random>

#include <Kernels/Utils.h>
#include <Kernels/KernelTraits.h>
#include <KdeTraits.h>
#include <Point2d.h>

namespace bbrcit {

// The one/two-point Epanechnikov product kernel in 2 dimensions 
// with bandwidth hx and hy is defined as follows:
//
// K(p;hx,hy) = C(hx,hy) * U(p,hx,hy)
// K(p,q;hx,hy) := C(hx,hy) * U(p,q,hx,hy)
//
// where: 
//
// + C(hx,hy) = (9 / 16) / (hx*hy)
//
// + U(p,hx,hy) = 
//       (1-px*px/(hx*hx)) * 
//         (1-py*py/(hy*hy)), if px*px/(hx*hx) < 1 and py*py/(hy*hy) < 1
//       0            , otherwise
//   
//
// + U(p,q,hx,hy) := U(p-q,hx,hy)
//
template<typename T=double>
class CUDA_ALIGN16 EpanechnikovProductKernel2d {

  public:

    using FloatType = T;

  public:
    CUDA_CALLABLE EpanechnikovProductKernel2d();
    CUDA_CALLABLE EpanechnikovProductKernel2d(T hx, T hy);

    CUDA_CALLABLE ~EpanechnikovProductKernel2d() = default;
    CUDA_CALLABLE EpanechnikovProductKernel2d(const EpanechnikovProductKernel2d<T>&) = default;
    CUDA_CALLABLE EpanechnikovProductKernel2d(EpanechnikovProductKernel2d<T>&&) = default;
    CUDA_CALLABLE EpanechnikovProductKernel2d& operator=(const EpanechnikovProductKernel2d<T>&) = default;
    CUDA_CALLABLE EpanechnikovProductKernel2d& operator=(EpanechnikovProductKernel2d<T>&&) = default;

    // evaluate C(hx,hy)
    CUDA_CALLABLE T normalization() const;

    // evaluate U(p,q,hx,hy)
    template<typename PointT>
      CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&) const;
    CUDA_CALLABLE T unnormalized_eval(const Point2d<T>&, const Point2d<T>&) const;

    // evaluate U(p,q,hx*a,hy*a)
    template<typename PointT>
      CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&, T) const;
    CUDA_CALLABLE T unnormalized_eval(const Point2d<T>&, const Point2d<T>&, T) const;

    // simulate a point from the kernel with local bandwidth correction `a`.
    // `e` is a random number engine from `std::random`. 
    // the result is stored in `p` with `p[i]` corresponding to component `i`. 
    template<typename RNG> 
      void simulate(RNG &e, std::vector<T> &p, T a = ConstantTraits<T>::one()) const;

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
EpanechnikovProductKernel2d<T>::EpanechnikovProductKernel2d() : hx_(1.0), hy_(1.0) {}

template<typename T>
EpanechnikovProductKernel2d<T>::EpanechnikovProductKernel2d(T hx, T hy) : hx_(hx), hy_(hy) {}

template<typename T>
inline T EpanechnikovProductKernel2d<T>::hx() const { return hx_; }

template<typename T>
inline T EpanechnikovProductKernel2d<T>::hy() const { return hy_; }

template<typename T>
inline void EpanechnikovProductKernel2d<T>::set_hx(T hx) { hx_ = hx; }

template<typename T>
inline void EpanechnikovProductKernel2d<T>::set_hy(T hy) { hy_ = hy; }

template<typename T>
inline void EpanechnikovProductKernel2d<T>::set_bandwidths(T hx, T hy) {
  hx_ = hx; hy_ = hy;
}

template <typename T>
inline T EpanechnikovProductKernel2d<T>::normalization() const {
  return EpanechnikovProduct2dTraits<T>::normalization(hx_, hy_);
}

template<typename T>
  template<typename PointT>
inline T EpanechnikovProductKernel2d<T>::unnormalized_eval(
    const PointT &p, const PointT &q) const {
  return unnormalized_eval(p, q, ConstantTraits<T>::one());
}

template<typename T>
inline T EpanechnikovProductKernel2d<T>::unnormalized_eval(
    const Point2d<T> &p, const Point2d<T> &q) const {
  return unnormalized_eval(p, q, ConstantTraits<T>::one());
}

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<typename T>
  template<typename PointT>
inline T EpanechnikovProductKernel2d<T>::unnormalized_eval(
    const PointT &lhs, const PointT &rhs, T a) const {
  return EpanechnikovProduct2dTraits<T>::kernel(
           (lhs[0]-rhs[0])*(lhs[0]-rhs[0])/(hx_*hx_*a*a)) *
         EpanechnikovProduct2dTraits<T>::kernel(
           (lhs[1]-rhs[1])*(lhs[1]-rhs[1])/(hy_*hy_*a*a));
}

template<typename T>
inline T EpanechnikovProductKernel2d<T>::unnormalized_eval(
    const Point2d<T> &lhs, const Point2d<T> &rhs, T a) const {
  return EpanechnikovProduct2dTraits<T>::kernel(
           (lhs.x()-rhs.x())*(lhs.x()-rhs.x())/(hx_*hx_*a*a)) *
         EpanechnikovProduct2dTraits<T>::kernel(
           (lhs.y()-rhs.y())*(lhs.y()-rhs.y())/(hy_*hy_*a*a));
}

template<typename T>
  template<typename RNG> 
void EpanechnikovProductKernel2d<T>::simulate(RNG &e, std::vector<T> &p, T a) const {

  static std::uniform_real_distribution<T> d(-1, 1);

  p.resize(2);
  p[0] = a * hx_ * epanechnikov_choice(d(e), d(e), d(e));
  p[1] = a * hy_ * epanechnikov_choice(d(e), d(e), d(e));
}

}

#endif
