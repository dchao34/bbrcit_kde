#ifndef BBRCITKDE_GAUSSIANKERNEL_H__
#define  BBRCITKDE_GAUSSIANKERNEL_H__

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>

#include <GeomUtils.h>

namespace bbrcit {

// The Gaussian kernel in D dimensions is the following:
//
// K(x) = unit_volume * exp{ -0.5 * x'x }
//
// where: 
//
// + x'x: 
//     dot product of x. 
// + unit_volume: 
//     (2 \pi)^{-D/2}
//
template<int D, typename T=double>
class GaussianKernel {

  public:

    static constexpr int dim() { return D; }

  public:
    GaussianKernel();
    GaussianKernel(T bandwidth);

    ~GaussianKernel() = default;
    GaussianKernel(const GaussianKernel&) = default;
    GaussianKernel(GaussianKernel&&) = default;
    GaussianKernel& operator=(const GaussianKernel&) = default;
    GaussianKernel& operator=(GaussianKernel&&) = default;

    // evaluate the value of the kernel at point p
    template<typename PointT>
    T eval(const PointT &p) const;

    // compute the normalization
    T normalization() const;

    // evaluate the kernel at point p, but do not include the 
    // normalization factor. e.g. evalutes exp{-0.5 * x'x}
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
GaussianKernel<D, T>::GaussianKernel() : bandwidth_(1.0) {}

template<int D, typename T>
GaussianKernel<D,T>::GaussianKernel(T bw) : bandwidth_(bw) {}

template<int D, typename T>
inline T GaussianKernel<D,T>::bandwidth() const { return bandwidth_; }

template<int D, typename T>
inline void GaussianKernel<D,T>::set_bandwidth(T bw) { bandwidth_ = bw; }

template<int D, typename T>
  template<typename PointT>
inline T GaussianKernel<D,T>::eval(const PointT &p) const {
  return normalization() * unnormalized_eval(p);
}

template<int D, typename T>
  template<typename PointT>
inline T GaussianKernel<D,T>::unnormalized_eval(const PointT &p) const {
  return exp(-0.5 * DotProduct(p, p) / (bandwidth_ * bandwidth_) );
}

template <int D, typename T>
inline T GaussianKernel<D,T>::normalization() const {
  return pow(2*M_PI, -D/2.0) / pow(bandwidth_, D);
}


// Specializations
// ---------------

// 1. normalization()
// ------------------

template <>
inline float GaussianKernel<1,float>::normalization() const {
  // 1 / (sqrt(2pi) * h); 1 / sqrt(2pi) = 0.39894228040143267794
  return 0.3989422804f / bandwidth_;
}

template <>
inline double GaussianKernel<1,double>::normalization() const {
  // 1 / (sqrt(2pi) * h); 1 / sqrt(2pi) = 0.39894228040143267794
  return 0.3989422804014326779 / bandwidth_;
}

template <>
inline float GaussianKernel<2,float>::normalization() const {
  // 1 / (2pi) * h * h; 1 / (2pi) = 0.15915494309189533577
  return 0.1591549431f / (bandwidth_ * bandwidth_);
}

template <>
inline double GaussianKernel<2,double>::normalization() const {
  // 1 / (2pi) * h * h; 1 / (2pi) = 0.15915494309189533577
  return 0.1591549430918953357 / (bandwidth_ * bandwidth_);
}

// 2. unnormalized_eval()
// ----------------------

template<>
  template<typename PointT>
inline float GaussianKernel<1,float>::unnormalized_eval(const PointT &p) const {
  return expf(-0.5f * DotProduct<PointT,float>(p, p) / (bandwidth_ * bandwidth_) );
}

template<>
  template<typename PointT>
inline float GaussianKernel<2,float>::unnormalized_eval(const PointT &p) const {
  return expf(-0.5f * DotProduct<PointT,float>(p, p) / (bandwidth_ * bandwidth_) );
}

template<>
  template<typename PointT>
inline double GaussianKernel<1,double>::unnormalized_eval(const PointT &p) const {
  return exp(-0.5 * DotProduct<PointT,double>(p, p) / (bandwidth_ * bandwidth_) );
}

template<>
  template<typename PointT>
inline double GaussianKernel<2,double>::unnormalized_eval(const PointT &p) const {
  return exp(-0.5 * DotProduct<PointT,double>(p, p) / (bandwidth_ * bandwidth_) );
}

}

#endif
