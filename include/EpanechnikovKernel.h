#ifndef BBRCITKDE_EPANECHNIKOVKERNEL_H__
#define  BBRCITKDE_EPANECHNIKOVKERNEL_H__

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>

#include <GeomUtils.h>

namespace bbrcit {

// The Epanechnikov kernel in D dimensions is the following:
//
// K(x) = 0.5 * unit_volume * (D+2) (1 - x'x)
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
    static const double normalization;

    // evaluate the value of the kernel at point p
    template<typename PointT>
    double eval(const PointT &p) const;

    // evaluate the kernel at point p, but do not include the 
    // normalization factor. e.g. evalutes (1-x'x), but does not 
    // multiply the leading constants
    template<typename PointT>
    double unnormalized_eval(const PointT&) const;

  private:
    static double compute_normalization();
    
};

// Implementations
// ---------------

template<int D>
  template<typename PointT>
double EpanechnikovKernel<D>::eval(const PointT &p) const {
  double dot_prod = DotProduct(p, p);
  return dot_prod < 1.0 ? (1-dot_prod) * normalization : 0; 
}

template<int D>
  template<typename PointT>
double EpanechnikovKernel<D>::unnormalized_eval(const PointT &p) const {
  double dot_prod = DotProduct(p, p);
  return dot_prod < 1.0 ? 1-dot_prod : 0;
}

template <int D>
const double EpanechnikovKernel<D>::normalization = EpanechnikovKernel<D>::compute_normalization();

template <int D>
inline double EpanechnikovKernel<D>::compute_normalization() {
  return 0.5 * (D+2) / (std::pow(M_PI, D/2.0) / std::tgamma(1+D/2.0));
}

}

#endif
