#ifndef BBRCITKDE_EPANECHNIKOVKERNEL_H__
#define BBRCITKDE_EPANECHNIKOVKERNEL_H__

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>

#include <KdeTraits.h>

namespace bbrcit {

// The Epanechnikov kernel in D dimensions is the following:
//
// K(x) = 
//     0.5 * unit_volume/h^D * (D+2) (1 - x'x/(h*h)), if x'x < 1
//     0                                              otherwise
//
// where: 
//
// + x'x: 
//     dot product of x. 
// + unit_volume: 
//     volume of a D-dim unit sphere, which is 
//     \pi^{D/2} * \Gamma(1+D/2), 
//     where Gamma is the gamma function. 
//
template<int D, typename T=double>
class EpanechnikovKernel {

  public:

    static constexpr int dim() { return D; }

  public:

    EpanechnikovKernel();
    EpanechnikovKernel(T bandwidth);

    ~EpanechnikovKernel() = default;
    EpanechnikovKernel(const EpanechnikovKernel&) = default;
    EpanechnikovKernel(EpanechnikovKernel&&) = default;
    EpanechnikovKernel& operator=(const EpanechnikovKernel&) = default;
    EpanechnikovKernel& operator=(EpanechnikovKernel&&) = default;

    // compute the normalization
    T normalization() const;

    // evaluate the two point kernel, but do not include the 
    // normalization factor. 
    template<typename PointT>
    T unnormalized_eval(const PointT&, const PointT&) const;

    // get/set the bandwidth
    T bandwidth() const;
    void set_bandwidth(T);

  private: 
    T bandwidth_;


    template<typename PointT>
    T point_arg_eval(const PointT&, const PointT&) const;
    
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

template<int D, typename T>
  template<typename PointT>
inline T EpanechnikovKernel<D,T>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return std::max(1.0 - point_arg_eval(p, q), 0.0);
}

template <int D, typename T>
inline T EpanechnikovKernel<D,T>::normalization() const {
  return 0.5 * (D+2) / (std::pow(M_PI, D/2.0) / std::tgamma(1+D/2.0)) 
                     / std::pow(bandwidth_, D);
}

// evaluates the (x-y)'(x-y)/h*h part for the kernel argument. 
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

// Specializations
// ---------------

// 1. normalization()
// ------------------

template <>
inline float EpanechnikovKernel<1,float>::normalization() const {
  // 0.5 * 3 / (sqrt(pi) * gamma(3/2)) / h = 0.75 / h
  return 0.75f / bandwidth_;
}

template <>
inline double EpanechnikovKernel<1,double>::normalization() const {
  // 0.5 * 3 / (sqrt(pi) * gamma(3/2)) / h = 0.75 / h
  return 0.75 / bandwidth_;
}

template <>
inline float EpanechnikovKernel<2,float>::normalization() const {
  // 0.5 * 4 / (pi*h*h) = 2 / (pi*h*h) =  0.63661977236758134307553 / (h*h)
  return 0.6366197723f / (bandwidth_ * bandwidth_);
}

template <>
inline double EpanechnikovKernel<2,double>::normalization() const {
  // 0.5 * 4 / (pi*h*h) = 2 / (pi*h*h) =  0.63661977236758134307553 / (h*h)
  return 0.6366197723675813431 / (bandwidth_ * bandwidth_);
}

// 2. unnormalized_eval()
// ----------------------

template<>
  template<typename PointT>
inline float EpanechnikovKernel<1,float>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return fmaxf(1.0f - point_arg_eval(p, q), 0.0f);
}

template<>
  template<typename PointT>
inline double EpanechnikovKernel<1,double>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return fmax(1.0 - point_arg_eval(p, q), 0.0);
}

template<>
  template<typename PointT>
inline float EpanechnikovKernel<2,float>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return fmaxf(1.0f - point_arg_eval(p, q), 0.0f);
}

template<>
  template<typename PointT>
inline double EpanechnikovKernel<2,double>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return fmax(1.0 - point_arg_eval(p, q), 0.0);
}


}

#endif
