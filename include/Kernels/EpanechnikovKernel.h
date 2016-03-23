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
#include <Kernels/Utils.h>
#include <KdeTraits.h>
#include <Point2d.h>
#include <Point1d.h>

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
//       1 - x'x/(h*h), if x'x/(h*h) < 1
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

    // evaluate U(x,y,h*a)
    template<typename PointT>
    CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&, T a) const;

    // simulate a point from the kernel with local bandwidth correction `a`.
    // `e` is a random number engine from `std::random`. 
    // the result is stored in `p` with `p[i]` corresponding to component `i`. 
    template<typename RNG> 
      void simulate(RNG &e, std::vector<T> &p, T a = ConstantTraits<T>::one()) const;

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
    CUDA_CALLABLE T point_arg_eval(const Point1d<T>&, const Point1d<T>&, T a) const;

    
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
  return unnormalized_eval(p, q, ConstantTraits<T>::one());
}

#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
template<int D, typename T>
  template<typename PointT>
inline T EpanechnikovKernel<D,T>::unnormalized_eval(const PointT &p, const PointT &q, T a) const {
  return EpanechnikovTraits<D,T>::kernel(point_arg_eval(p, q, a));
}

template <int D, typename T>
inline T EpanechnikovKernel<D,T>::normalization() const {
  return EpanechnikovTraits<D,T>::normalization(bandwidth_);
}

template<int D, typename T>
  template<typename PointT>
T EpanechnikovKernel<D,T>::point_arg_eval(const PointT &lhs, const PointT &rhs, T a) const {

  T result = ConstantTraits<T>::zero(); 

  T diff; 
  for (int i = 0; i < D; ++i) {
    diff = lhs[i] - rhs[i];
    result += diff * diff;
  }

  return result / (bandwidth_ * bandwidth_ * a * a);
}

template<int D, typename T>
inline T EpanechnikovKernel<D,T>::point_arg_eval(
    const Point2d<T> &lhs, const Point2d<T> &rhs, T a) const {
  return ((lhs.x()-rhs.x())*(lhs.x()-rhs.x()) +
          (lhs.y()-rhs.y())*(lhs.y()-rhs.y())) / (bandwidth_*bandwidth_*a*a);
}

template<int D, typename T>
inline T EpanechnikovKernel<D,T>::point_arg_eval(
    const Point1d<T> &lhs, const Point1d<T> &rhs, T a) const {
  return ((lhs.x()-rhs.x())*(lhs.x()-rhs.x())) / (bandwidth_*bandwidth_*a*a);
}


template<>
  template<typename RNG> 
void EpanechnikovKernel<1, float>::simulate(RNG &e, std::vector<float> &p, float a) const {

  static std::uniform_real_distribution<float> d(-1, 1);
  p.resize(1);
  p[0] = a * bandwidth_ * epanechnikov_choice(d(e), d(e), d(e));

}

template<>
  template<typename RNG> 
void EpanechnikovKernel<1, double>::simulate(RNG &e, std::vector<double> &p, double a) const {

  static std::uniform_real_distribution<double> d(-1, 1);
  p.resize(1);
  p[0] = a * bandwidth_ * epanechnikov_choice(d(e), d(e), d(e));

}


}

#endif
