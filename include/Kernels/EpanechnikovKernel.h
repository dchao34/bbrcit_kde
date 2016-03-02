#ifndef BBRCITKDE_EPANECHNIKOVKERNEL_H__
#define BBRCITKDE_EPANECHNIKOVKERNEL_H__

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>

#include <GeomUtils.h>

namespace bbrcit {

// The Epanechnikov kernel in D dimensions is the following:
//
// K(x) = 
//     0.5 * unit_volume * (D+2) (1 - x'x), if x'x < 1
//     0                                    otherwise
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

    // evaluate the value of the kernel at point p
    template<typename PointT>
    T eval(const PointT &p) const;

    // compute the normalization
    T normalization() const;

    // evaluate the kernel at point p, but do not include the 
    // normalization factor. e.g. evalutes (1-x'x), but does not 
    // multiply the leading constants
    template<typename PointT>
    T unnormalized_eval(const PointT&) const;

    // get/set the bandwidth
    T bandwidth() const;
    void set_bandwidth(T);

  private: 
    T bandwidth_;

    
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
inline T EpanechnikovKernel<D,T>::eval(const PointT &p) const {
  return normalization() * unnormalized_eval(p); 
}

template<int D, typename T>
  template<typename PointT>
inline T EpanechnikovKernel<D,T>::unnormalized_eval(const PointT &p) const {
  return std::max(1.0 - DotProduct(p,p) / (bandwidth_ * bandwidth_), 0.0);
}

template <int D, typename T>
inline T EpanechnikovKernel<D,T>::normalization() const {
  return 0.5 * (D+2) / (std::pow(M_PI, D/2.0) / std::tgamma(1+D/2.0)) 
                     / std::pow(bandwidth_, D);
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
inline float EpanechnikovKernel<1,float>::unnormalized_eval(const PointT &p) const {
  return fmaxf(1.0f - DotProduct<PointT,float>(p,p)/(bandwidth_*bandwidth_), 0.0f);
}

template<>
  template<typename PointT>
inline float EpanechnikovKernel<2,float>::unnormalized_eval(const PointT &p) const {
  return fmaxf(1.0f - DotProduct<PointT,float>(p,p)/(bandwidth_*bandwidth_), 0.0f);
}

template<>
  template<typename PointT>
inline double EpanechnikovKernel<1,double>::unnormalized_eval(const PointT &p) const {
  return fmax(1.0 - DotProduct<PointT,double>(p,p)/(bandwidth_*bandwidth_), 0.0);
}

template<>
  template<typename PointT>
inline double EpanechnikovKernel<2,double>::unnormalized_eval(const PointT &p) const {
  return fmax(1.0 - DotProduct<PointT,double>(p,p)/(bandwidth_*bandwidth_), 0.0);
}


}

#endif
