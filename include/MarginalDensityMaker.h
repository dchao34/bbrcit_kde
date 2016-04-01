#ifndef BBRCIT_MARGINALDENSITYMAKER_H__
#define BBRCIT_MARGINALDENSITYMAKER_H__

#include <vector>

#include <Kernels/GaussianProductKernel2d.h>
#include <Kernels/GaussianKernel.h>
#include <Kernels/EpanechnikovProductKernel2d.h>
#include <Kernels/EpanechnikovKernel.h>
#include <KernelDensity.h>
#include <Attributes/AdaKdeAttributes.h>

namespace bbrcit {

// Point Marginalizers
// -------------------

// Input:
//   + `pt2d`: 2d DecoratedPoint whose attributes has at least all 
//              those in AdaKdeAttributes.
//   + `dim_idx`: dimension index to marginalize. 0: x, 1: y.
template <typename FloatT, typename AttrT>
DecoratedPoint<1,AttrT,FloatT> marginalize_decorated_point2d(
    const DecoratedPoint<2,AttrT,FloatT> &pt2d, 
    size_t dim_idx) {

  DecoratedPoint<1,AttrT,FloatT> pt1d;

  // marginalize the point coordinate
  pt1d[0] = pt2d[dim_idx];

  // transcribe the attributes. mass is just weight / abw with no raised powers. 
  pt1d.attributes().set_weight(pt2d.attributes().weight());
  pt1d.attributes().set_lower_abw(pt2d.attributes().lower_abw());
  pt1d.attributes().set_upper_abw(pt2d.attributes().upper_abw());
  pt1d.attributes().set_mass(
      pt1d.attributes().weight() / pt1d.attributes().abw()
  );
  return pt1d;
}


// Kernel Marginalizers
// --------------------

// class with metafunctions that define how kernels map to their marginal kernels
template <typename T>
class MarginalKernelAssociator {};

// partial specialization for EpanechnikovProductKernel2d
template <typename T>
class MarginalKernelAssociator<EpanechnikovProductKernel2d<T>> {

  public:
    using KernelType = EpanechnikovProductKernel2d<T>;
    using Kernel1dType = EpanechnikovKernel<1,T>;

    static Kernel1dType make_marginal_kernel(const KernelType &kern, size_t dim_idx) {
      Kernel1dType kern1d(dim_idx == 0 ? kern.hx() : kern.hy());
      return kern1d;
    }
};

// partial specialization for GaussianProductKernel2d
template <typename T>
class MarginalKernelAssociator<GaussianProductKernel2d<T>> {

  public:
    using KernelType = GaussianProductKernel2d<T>;
    using Kernel1dType = GaussianKernel<1,T>;

    static Kernel1dType make_marginal_kernel(const KernelType &kern, size_t dim_idx) {
      Kernel1dType kern1d(dim_idx == 0 ? kern.hx() : kern.hy());
      return kern1d;
    }
};


// Density Marginalizers
// ---------------------

// given a kernel density in 2d, return the corresponding marginal density in 1d. 
template<typename KernelT, typename FloatT, typename AttrT>
KernelDensity<1,typename MarginalKernelAssociator<KernelT>::Kernel1dType,FloatT,AttrT>
make_marginal_density(const KernelDensity<2,KernelT,FloatT,AttrT> &kde, size_t dim_idx) {

  using Kernel1dType = typename MarginalKernelAssociator<KernelT>::Kernel1dType;
  using Kde1dType = KernelDensity<1,Kernel1dType,FloatT,AttrT>;
  using Point1dType = typename Kde1dType::DataPointType;

  // marginalize the data points
  std::vector<Point1dType> points1d;
  for (const auto &p : kde.points()) {
    points1d.push_back(marginalize_decorated_point2d(p, dim_idx));
  }

  // construct the corresponding 1d density estimator
  Kde1dType kde1d(std::move(points1d), kde.data_tree().leaf_nmax());

  // configure the 1d kernel 
  kde1d.kernel() = MarginalKernelAssociator<KernelT>::make_marginal_kernel(kde.kernel(), dim_idx);

  return kde1d;
}

// convenience wrapper for `make_marginal_density`
template<typename KernelT, typename FloatT, typename AttrT>
inline KernelDensity<1,typename MarginalKernelAssociator<KernelT>::Kernel1dType,FloatT,AttrT>
make_marginal_density_x(const KernelDensity<2,KernelT,FloatT,AttrT> &kde) {
  return make_marginal_density(kde, 0);
}

// convenience wrapper for `make_marginal_density`
template<typename KernelT, typename FloatT, typename AttrT>
inline KernelDensity<1,typename MarginalKernelAssociator<KernelT>::Kernel1dType,FloatT,AttrT>
make_marginal_density_y(const KernelDensity<2,KernelT,FloatT,AttrT> &kde) {
  return make_marginal_density(kde, 1);
}


}

#endif
