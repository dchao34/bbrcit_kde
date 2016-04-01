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
    using MarginalKernelType = EpanechnikovKernel<1,T>;

    static MarginalKernelType marginal_kernel(const KernelType &kern, size_t dim_idx) {
      MarginalKernelType marginal_kern(dim_idx == 0 ? kern.hx() : kern.hy());
      return marginal_kern;
    }
};

// partial specialization for GaussianProductKernel2d
template <typename T>
class MarginalKernelAssociator<GaussianProductKernel2d<T>> {

  public:
    using KernelType = GaussianProductKernel2d<T>;
    using MarginalKernelType = GaussianKernel<1,T>;

    static MarginalKernelType marginal_kernel(const KernelType &kern, size_t dim_idx) {
      MarginalKernelType marginal_kern(dim_idx == 0 ? kern.hx() : kern.hy());
      return marginal_kern;
    }
};


// Density Marginalizers
// ---------------------

// class with metafunctions that define how densities kernels map to their marginal densities
template <typename KdeT, 
          typename KernelT = typename KdeT::KernelType, 
          typename FloatT= typename KdeT::FloatType, 
          typename AttrT= typename KdeT::AttributeType>
class MarginalDensityAssociator {};

// specialization for KernelDensity<2,...>
template <typename KernelT, typename FloatT, typename AttrT>
class MarginalDensityAssociator<KernelDensity<2,KernelT,FloatT,AttrT>, 
                                KernelT,
                                FloatT, AttrT> {
  public:
    using KernelDensityType = KernelDensity<2,KernelT,FloatT,AttrT>;
    using MarginalKernelType = typename MarginalKernelAssociator<KernelT>::MarginalKernelType;
    using MarginalDensityType = KernelDensity<1,MarginalKernelType,FloatT,AttrT>;
    using Point1d = typename MarginalDensityType::DataPointType;
    
    // given a kernel density in 2d, return the corresponding marginal density in 1d. 
    static MarginalDensityType marginal_density(const KernelDensityType &kde, size_t dim_idx) {

      // marginalize the data points
      std::vector<Point1d> points;
      for (const auto &p : kde.points()) {
        points.push_back(marginalize_decorated_point2d(p, dim_idx));
      }

      // construct the corresponding 1d density estimator
      MarginalDensityType marginal_kde(std::move(points), kde.data_tree().leaf_nmax());

      // configure the marginal kernel 
      marginal_kde.kernel() = MarginalKernelAssociator<KernelT>::marginal_kernel(kde.kernel(), dim_idx);

      return marginal_kde;
    }

    // convenience wrapper for making x marginals
    static MarginalDensityType marginal_density_x(const KernelDensityType &kde) {
      return marginal_density(kde, 0);
    }

    // convenience wrapper for making y marginals
    static MarginalDensityType marginal_density_y(const KernelDensityType &kde) {
      return marginal_density(kde, 1);
    }

};

}

#endif
