#ifndef BBRCIT_MARGINALDENSITYMAKER_H__
#define BBRCIT_MARGINALDENSITYMAKER_H__

#include <vector>

#include <Kernels/EpanechnikovProductKernel2d.h>
#include <KernelDensity.h>
#include <Attributes/AdaKdeAttributes.h>

namespace bbrcit {

template <typename Point1dT, typename Point2dT>
Point1dT project_point1d(const Point2dT &p2d, size_t dim_idx) {
  Point1dT p1d;
  p1d[0] = p2d[dim_idx];
  p1d.attributes().set_weight(p2d.attributes().weight());
  p1d.attributes().set_lower_abw(p2d.attributes().lower_abw());
  p1d.attributes().set_upper_abw(p2d.attributes().upper_abw());
  p1d.attributes().set_mass(
      p1d.attributes().weight() / p1d.attributes().abw()
  );
  return p1d;
}


template<typename KFloatT, typename FloatT>
class EpanechnikovProduct2dKdeMarginalizer {

  public:

    using FloatType = FloatT;
    using KFloatType = KFloatT;
    using AttributeType = AdaKdeAttributes<FloatT>;

    using KernelType = EpanechnikovProductKernel2d<KFloatT>;
    using KernelDensityType = KernelDensity<2, KernelType, FloatT, AttributeType>;
    using DataPointType = typename KernelDensityType::DataPointType;

    using Kernel1dType = EpanechnikovKernel<1,KFloatT>;
    using KernelDensity1dType = KernelDensity<1, Kernel1dType, FloatT, AttributeType>;
    using DataPoint1dType = typename KernelDensity1dType::DataPointType;

    static KFloatType get_bandwidth(const KernelDensityType &kde, size_t dim_idx) {
      return dim_idx == 0 ? kde.kernel().hx() : kde.kernel().hy();
    }

    static KernelDensity1dType make_1d_marginal_density(
        const KernelDensityType &kde, size_t dim_idx, int leaf_nmax) {

      std::vector<DataPoint1dType> point1d;
      for (const auto &p : kde.points()) {
        point1d.push_back(project_point1d<DataPoint1dType>(p, dim_idx));
      }
      KernelDensity1dType kde1d(std::move(point1d), leaf_nmax);

      kde1d.kernel().set_bandwidth(get_bandwidth(kde, dim_idx));

      return kde1d;
    }

    static KernelDensity1dType make_marginal_density_x(
        const KernelDensityType &kde) {
      return make_1d_marginal_density(kde, 0, kde.data_tree().leaf_nmax());
    }

    static KernelDensity1dType make_marginal_density_y(
        const KernelDensityType &kde) {
      return make_1d_marginal_density(kde, 1, kde.data_tree().leaf_nmax());
    }

};



}

#endif
