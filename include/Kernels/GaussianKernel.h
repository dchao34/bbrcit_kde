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
template<int D>
class GaussianKernel {

  public:

    static constexpr int dim() { return D; }

  public:
    GaussianKernel();
    GaussianKernel(double bandwidth);

    ~GaussianKernel() = default;
    GaussianKernel(const GaussianKernel&) = default;
    GaussianKernel(GaussianKernel&&) = default;
    GaussianKernel& operator=(const GaussianKernel&) = default;
    GaussianKernel& operator=(GaussianKernel&&) = default;

    // evaluate the value of the kernel at point p
    template<typename PointT>
    double eval(const PointT &p) const;

    // compute the normalization
    double normalization() const;

    // evaluate the kernel at point p, but do not include the 
    // normalization factor. e.g. evalutes exp{-0.5 * x'x}
    template<typename PointT>
    double unnormalized_eval(const PointT&) const;

    // get/set the bandwidth
    double bandwidth() const;
    void set_bandwidth(double);

  private:
    double bandwidth_;
    
};

// Implementations
// ---------------

template<int D>
GaussianKernel<D>::GaussianKernel() : bandwidth_(1.0) {}

template<int D>
GaussianKernel<D>::GaussianKernel(double bw) : bandwidth_(bw) {}

template<int D>
inline double GaussianKernel<D>::bandwidth() const { return bandwidth_; }

template<int D>
inline void GaussianKernel<D>::set_bandwidth(double bw) { bandwidth_ = bw; }

template<int D>
  template<typename PointT>
inline double GaussianKernel<D>::eval(const PointT &p) const {
  return normalization() * unnormalized_eval(p);
}

template<int D>
  template<typename PointT>
inline double GaussianKernel<D>::unnormalized_eval(const PointT &p) const {
  return std::exp(-0.5 * DotProduct(p, p) / (bandwidth_ * bandwidth_) );
}

template <int D>
inline double GaussianKernel<D>::normalization() const {
  return std::pow(2*M_PI, -D/2.0) / std::pow(bandwidth_, D);
}

}

#endif
