#ifndef BBRCITKDE_EPANECHNIKOVPRODUCTKERNEL_H__
#define BBRCITKDE_EPANECHNIKOVPRODUCTKERNEL_H__

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <initializer_list>
#include <vector>

#include <GeomUtils.h>

namespace bbrcit {

// The Epanechnikov product kernel in D dimensions is the following:
//
// K(x) = 
//     \Prod^{D}_{i=1} 0.75/h_i * (1 - x_i^2 / h_i^2), if x_i^2/h_i^2 < 1 for all i, 
//     0                                           otherwise
//
template<int D>
class EpanechnikovProductKernel {

  public:

    static constexpr int dim() { return D; }

  public:

    EpanechnikovProductKernel();
    EpanechnikovProductKernel(const std::initializer_list<double>&);
    ~EpanechnikovProductKernel() = default;
    EpanechnikovProductKernel(const EpanechnikovProductKernel&) = default;
    EpanechnikovProductKernel(EpanechnikovProductKernel&&) = default;
    EpanechnikovProductKernel& operator=(const EpanechnikovProductKernel&) = default;
    EpanechnikovProductKernel& operator=(EpanechnikovProductKernel&&) = default;

    // evaluate the value of the kernel at point p
    template<typename PointT>
    double eval(const PointT &p) const;

    // compute the normalization
    double normalization() const;

    // evaluate the kernel at point p, but do not include the 
    // normalization factor. e.g. evalutes exp{-0.5 * x'x}
    template<typename PointT>
    double unnormalized_eval(const PointT&) const;

    // index access to bandwidths
    double bandwidth(size_t) const;

    // set bandwidths
    void set_bandwidth(size_t, double);
    void set_bandwidths(const std::initializer_list<double> &);

  private:
    std::vector<double> bandwidths_;
    
};

// Implementations
// ---------------

template<int D>
EpanechnikovProductKernel<D>::EpanechnikovProductKernel() : bandwidths_(D, 1.0) {}

template<int D>
EpanechnikovProductKernel<D>::EpanechnikovProductKernel(
    const std::initializer_list<double> &li): bandwidths_(li) {
  if (bandwidths_.size() != D) { bandwidths_.resize(D); }
}

template<int D>
inline double EpanechnikovProductKernel<D>::bandwidth(size_t i) const { 
  return bandwidths_[i]; 
}

template<int D>
inline void EpanechnikovProductKernel<D>::set_bandwidth(size_t i, double bw) { 
  bandwidths_[i] = bw; 
}

template<int D>
void EpanechnikovProductKernel<D>::set_bandwidths(
    const std::initializer_list<double> &li) {
  bandwidths_ = li;
  if (bandwidths_.size() != D) { bandwidths_.resize(D); }
}

template<int D>
  template<typename PointT>
inline double EpanechnikovProductKernel<D>::eval(const PointT &p) const {
  return normalization() * unnormalized_eval(p);
}

template<int D>
  template<typename PointT>
double EpanechnikovProductKernel<D>::unnormalized_eval(const PointT &p) const {

  double result = 1.0;

  double dot_prod;
  for (int i = 0; i < D; ++i) {
    dot_prod = p[i] * p[i] / (bandwidths_[i] * bandwidths_[i]);
    dot_prod < 1 ? result *= (1 - dot_prod) : result = 0;
  }

  return result;
}

template <int D>
double EpanechnikovProductKernel<D>::normalization() const {
  double result = std::pow(0.75, D);
  for (const auto bw : bandwidths_) { result /= bw; }
  return result;
}

}

#endif
