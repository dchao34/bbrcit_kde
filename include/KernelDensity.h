#ifndef BBRCITKDE_KERNELDENSITY_H__
#define BBRCITKDE_KERNELDENSITY_H__

#include <limits>
#include <queue>
#include <stack>
#include <tuple>

#include <Kdtree.h>
#include <DecoratedPoint.h>
#include <PointWeights.h>
#include <EpanechnikovKernel.h>
#include <KdeTraits.h>

namespace bbrcit {

template<int D, typename KT, typename PT, typename NT, typename FT> class KernelDensity;

template<int D, typename KT, typename PT, typename NT, typename FT>
void swap(KernelDensity<D,KT,PT,NT,FT>&, KernelDensity<D,KT,PT,NT,FT>&);


// KernelDensity<> implements a kernel density estimator in D dimensions using 
// a D dimensional Kdtree. 
template<int D, 
         typename KernelT=EpanechnikovKernel<D>,
         typename PAttrT=PointWeights<double>,
         typename NAttrT=PAttrT,
         typename FloatT=double>
class KernelDensity : public Kdtree<D,PAttrT,NAttrT,FloatT> {

  private:
    using KdtreeType = Kdtree<D,PAttrT,NAttrT,FloatT>;
    using KernelDensityType = KernelDensity<D,KernelT,PAttrT,NAttrT,FloatT>;
    using NodeType = typename KdtreeType::Node;

  public: 
    using KernelType = KernelT;
    using DataPointType = typename KdtreeType::DataPointType;
    using RectangleType = typename KdtreeType::RectangleType;
    using GeomPointType = typename DataPointType::PointType;
    using FloatType = FloatT;

    friend void swap<>(KernelDensityType&, KernelDensityType&);

    static constexpr int dim() { return D; }

    // default constructor: estimator on an empty set. 
    KernelDensity();

    // construct a kernel density estimator with bandwidth `bw` over the points `data`.
    KernelDensity(const std::vector<DataPointType> &data, FloatType bw, int leaf_nmax=2);
    KernelDensity(std::vector<DataPointType>&&, FloatType bw, int leaf_nmax=2);

    // copy-control
    KernelDensity(const KernelDensityType&) = default;
    KernelDensity(KernelDensityType&&) noexcept = default;
    KernelDensityType& operator=(const KernelDensityType&) = default;
    KernelDensityType& operator=(KernelDensityType&&) = default;
    ~KernelDensity();

    // get/set bandwidth
    FloatType bandwidth() const;
    void set_bandwidth(FloatType);

    // evaluate the kde at point `p`. relative error is no more than `rtol`.
    FloatT eval(const GeomPointType &p, FloatType rtol) const;

    // mainly for debugging: naive kde evaluation; slow... O(n^2). 
    FloatT naive_eval(const GeomPointType&) const;

  private:

    // bandwidth of the estimator
    FloatType bandwidth_;

    // the kernel function of the estimator 
    KernelType kernel_;

    void single_tree(const GeomPointType&, FloatType, 
                     NodeType*, const RectangleType&, size_t, 
                     FloatType&, FloatType&) const;

};

template<int D, typename KT, typename PT, typename NT, typename FT>
typename KernelDensity<D,KT,PT,NT,FT>::FloatType
KernelDensity<D,KT,PT,NT,FT>::eval(
    const GeomPointType &p, FloatType rtol) const {

  FloatType upper_bound = (this->root_)->attr_.get_weight(); 
  FloatType lower_bound = ConstantTraits<FloatType>::zero();

  this->single_tree(p, rtol, this->root_, this->bbox_, 0, upper_bound, lower_bound);

  FloatType result = lower_bound + (upper_bound - lower_bound) / 2;
  result *= KernelType::normalization;
  result /= (std::pow(bandwidth_, D) * this->points_.size());
  return result;
}



template<int D, typename KT, typename PT, typename NT, typename FT>
void KernelDensity<D,KT,PT,NT,FT>::single_tree(
    const GeomPointType &p, FloatType rtol, 
    NodeType *r, const RectangleType &bounding_box, size_t depth,
    FloatType &u, FloatType &l) const {

  if (r == nullptr) { return; }

  GeomPointType proxy; 

  proxy[0] = bounding_box.min_dist(p);
  FloatType max_val = kernel_.unnormalized_eval(proxy/bandwidth_);
  FloatType du = r->attr_.get_weight() * (max_val - 1);

  proxy[0] = bounding_box.max_dist(p);
  FloatType min_val = kernel_.unnormalized_eval(proxy/bandwidth_);
  FloatType dl = r->attr_.get_weight() * min_val;

  if (max_val < std::numeric_limits<FloatType>::epsilon() ||
      std::abs(u+du-l-dl)/std::abs(l+dl) < rtol) {
    u += du; l += dl; return;
  }

  if (r->is_leaf()) {
    for (auto i = r->start_idx_; i <= r->end_idx_; ++i) {
      FloatType delta = kernel_.unnormalized_eval(
          (p - this->points_[i].point()) / this->bandwidth_
      );
      u+=delta; l+= delta;
    }
    u -= r->attr_.get_weight();
  } else {
    auto lower_bbox = bounding_box.lower_halfspace(depth, r->split_);
    auto upper_bbox = bounding_box.upper_halfspace(depth, r->split_);
    proxy[0] = lower_bbox.min_dist(p);
    FloatType lower_dist = kernel_.unnormalized_eval(proxy/bandwidth_);
    proxy[0] = upper_bbox.min_dist(p);
    FloatType upper_dist = kernel_.unnormalized_eval(proxy/bandwidth_);
    if (lower_dist < upper_dist) {
      this->single_tree(p, rtol, r->left_, lower_bbox, (depth+1)%D, u, l);
      this->single_tree(p, rtol, r->right_, upper_bbox, (depth+1)%D, u, l);
    } else {
      this->single_tree(p, rtol, r->right_, upper_bbox, (depth+1)%D, u, l);
      this->single_tree(p, rtol, r->left_, lower_bbox, (depth+1)%D, u, l);
    }
  }

}

template<int D, typename KT, typename PT, typename NT, typename FT>
typename KernelDensity<D,KT,PT,NT,FT>::FloatType
KernelDensity<D,KT,PT,NT,FT>::naive_eval(const GeomPointType &p) const {
  FloatType total = ConstantTraits<FloatType>::zero();
  for (const auto &datum : this->points_) {
    total += kernel_.unnormalized_eval( (p - datum.point()) / bandwidth_  );
  }
  total *= KernelType::normalization;
  total /= (std::pow(bandwidth_, D) * this->points_.size());
  return total;
}

template<int D, typename KT, typename PT, typename NT, typename FT>
inline typename KernelDensity<D,KT,PT,NT,FT>::FloatType
KernelDensity<D,KT,PT,NT,FT>::bandwidth() const { return bandwidth_; }

template<int D, typename KT, typename PT, typename NT, typename FT>
inline void KernelDensity<D,KT,PT,NT,FT>::set_bandwidth(FloatType bw) { 
  bandwidth_ = bw; 
}

template<int D, typename KT, typename PT, typename NT, typename FT>
void swap(KernelDensity<D,KT,PT,NT,FT> &lhs, KernelDensity<D,KT,PT,NT,FT> &rhs) {
  using KdtreeType = typename KernelDensity<D,KT,PT,NT,FT>::KdtreeType;
  using std::swap;
  swap(static_cast<KdtreeType&>(lhs), static_cast<KdtreeType&>(rhs));
  swap(lhs.bandwidth_, rhs.bandwidth_);
  swap(lhs.kernel_, rhs.kernel_);
  return;
}

template<int D, typename KT, typename PT, typename NT, typename FT>
KernelDensity<D,KT,PT,NT,FT>::KernelDensity() : KdtreeType(), bandwidth_(1), kernel_() {}

template<int D, typename KT, typename PT, typename NT, typename FT>
KernelDensity<D,KT,PT,NT,FT>::KernelDensity(
    const std::vector<DataPointType> &pts, FloatType bw, int leaf_max) 
  : KdtreeType(pts, leaf_max), bandwidth_(bw), kernel_() {}

template<int D, typename KT, typename PT, typename NT, typename FT>
KernelDensity<D,KT,PT,NT,FT>::KernelDensity(
    std::vector<DataPointType> &&pts, FloatType bw, int leaf_max) 
  : KdtreeType(std::move(pts), leaf_max), bandwidth_(bw), kernel_() {}

template<int D, typename KT, typename PT, typename NT, typename FT>
KernelDensity<D,KT,PT,NT,FT>::~KernelDensity() {}

}

#endif
