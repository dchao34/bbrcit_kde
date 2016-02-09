#ifndef BBRCITKDE_GAUSSIANPRODUCTKERNEL_H__
#define  BBRCITKDE_GAUSSIANPRODUCTKERNEL_H__

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <initializer_list>
#include <vector>

#include <GeomUtils.h>

namespace bbrcit {

// The Gaussian product kernel in D dimensions is the following:
//
// K(x) = \Prod^{D}_{i=1} (2\pi * h_i)^{-0.5} exp{ -0.5 * x_i^2 / h_i^2 }
//
template<int D>
class GaussianProductKernel {

  public:

    static constexpr int dim() { return D; }

  public:

    GaussianProductKernel();
    GaussianProductKernel(const std::initializer_list<double>&);
    ~GaussianProductKernel() = default;
    GaussianProductKernel(const GaussianProductKernel&) = default;
    GaussianProductKernel(GaussianProductKernel&&) = default;
    GaussianProductKernel& operator=(const GaussianProductKernel&) = default;
    GaussianProductKernel& operator=(GaussianProductKernel&&) = default;

    // evaluate the value of the kernel at point p
    template<typename PointT>
    double eval(const PointT &p) const;

    // compute the normalization
    double normalization() const;

    // evaluate the kernel at point p, but do not include the 
    // normalization factor. e.g. evalutes exp{-0.5 * x'x}
    template<typename PointT>
    double unnormalized_eval(const PointT&) const;

  private:
    std::vector<double> bandwidths_;
    
};

// Implementations
// ---------------

template<int D>
GaussianProductKernel<D>::GaussianProductKernel() : bandwidths_(D, 1.0) {}

template<int D>
GaussianProductKernel<D>::GaussianProductKernel(
    const std::initializer_list<double> &li): bandwidths_(li) {
  if (bandwidths_.size() != D) { bandwidths_.resize(D); }
}

template<int D>
  template<typename PointT>
inline double GaussianProductKernel<D>::eval(const PointT &p) const {
  return normalization() * unnormalized_eval(p);
}

template<int D>
  template<typename PointT>
double GaussianProductKernel<D>::unnormalized_eval(const PointT &p) const {
  double arg = 0.0;
  for (int i = 0; i < D; ++i) {
    arg += p[i]*p[i] / (bandwidths_[i]*bandwidths_[i]); 
  }
  return std::exp(-0.5 * arg);
}

template <int D>
double GaussianProductKernel<D>::normalization() const {
  double result = std::pow(2*M_PI, -D/2.0);
  for (const auto bw : bandwidths_) { result /= bw; }
  return result;
}

}

#endif
