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
template<int D>
class EpanechnikovKernel {

  public:

    static constexpr int dim() { return D; }

  public:

    EpanechnikovKernel();
    EpanechnikovKernel(double bandwidth);

    ~EpanechnikovKernel() = default;
    EpanechnikovKernel(const EpanechnikovKernel&) = default;
    EpanechnikovKernel(EpanechnikovKernel&&) = default;
    EpanechnikovKernel& operator=(const EpanechnikovKernel&) = default;
    EpanechnikovKernel& operator=(EpanechnikovKernel&&) = default;

    // evaluate the value of the kernel at point p
    template<typename PointT>
    double eval(const PointT &p) const;

    // compute the normalization
    double normalization() const;

    // evaluate the kernel at point p, but do not include the 
    // normalization factor. e.g. evalutes (1-x'x), but does not 
    // multiply the leading constants
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
EpanechnikovKernel<D>::EpanechnikovKernel() : bandwidth_(1.0) {}

template<int D>
EpanechnikovKernel<D>::EpanechnikovKernel(double bw) : bandwidth_(bw) {}

template<int D>
inline double EpanechnikovKernel<D>::bandwidth() const { return bandwidth_; }

template<int D>
inline void EpanechnikovKernel<D>::set_bandwidth(double bw) { bandwidth_ = bw; }

template<int D>
  template<typename PointT>
inline double EpanechnikovKernel<D>::eval(const PointT &p) const {
  return normalization() * unnormalized_eval(p); 
}

template<int D>
  template<typename PointT>
double EpanechnikovKernel<D>::unnormalized_eval(const PointT &p) const {
  double dot_prod = DotProduct(p, p) / (bandwidth_ * bandwidth_);
  return dot_prod < 1.0 ? 1-dot_prod : 0;
}

template <int D>
inline double EpanechnikovKernel<D>::normalization() const {
  return 0.5 * (D+2) / (std::pow(M_PI, D/2.0) / std::tgamma(1+D/2.0)) 
                     / std::pow(bandwidth_, D);
}

}

#endif
