#ifndef BBRCITKDE_KERNELDENSITY_H__
#define BBRCITKDE_KERNELDENSITY_H__

#include <limits>
#include <queue>
#include <stack>
#include <tuple>
#include <utility>

#include <Kdtree.h>
#include <DecoratedPoint.h>
#include <PointWeights.h>
#include <QueryTreeAttributes.h>
#include <EpanechnikovKernel.h>
#include <KdeTraits.h>

namespace bbrcit {

template<int D, typename KT, typename DT, typename QT, typename FT> class KernelDensity;

template<int D, typename KT, typename DT, typename QT, typename FT>
void swap(KernelDensity<D,KT,DT,QT,FT>&, KernelDensity<D,KT,DT,QT,FT>&);


// KernelDensity<> implements a kernel density estimator in D dimensions using 
// a D dimensional Kdtree. 
template<int D, 
         typename KernelT=EpanechnikovKernel<D>,
         typename DAttrT=PointWeights<double>,
         typename QAttrT=QueryTreeAttributes<double>,
         typename FloatT=double>
class KernelDensity {

  private:

    using KernelDensityType = KernelDensity<D,KernelT,DAttrT,QAttrT,FloatT>;

    using DataTreeType = Kdtree<D,DAttrT,DAttrT,FloatT>; 
    using DataNodeType = typename DataTreeType::Node;

  public: 
    using DataPointType = typename DataTreeType::DataPointType;
    using GeomRectangleType = typename DataTreeType::RectangleType;
    using GeomPointType = typename DataPointType::PointType;

    using KernelType = KernelT;
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


    // primarily for debugging:
    bool empty() const { return data_tree_.root_ == nullptr; }
    size_t size() const { return data_tree_.size(); }
    const std::vector<DataPointType>& data_points() const { return data_tree_.points(); }
    void report_leaves(std::vector<std::pair<size_t,size_t>> &l) const { data_tree_.report_leaves(l); }
    void range_search(
        const GeomRectangleType &q, 
        std::vector<DataPointType> &r) const { data_tree_.range_search(q, r); }
    const DAttrT& root_attributes() const { return data_tree_.root_->attr_; }

    // mainly for debugging: naive kde evaluation; slow... O(n^2). 
    FloatT naive_eval(const GeomPointType&) const;

  private:

    // bandwidth of the estimator
    FloatType bandwidth_;

    // the kernel function of the estimator 
    KernelType kernel_;

    // the kernel function of the estimator 
    DataTreeType data_tree_;

    void single_tree(const GeomPointType&, FloatType, 
                     DataNodeType*, const GeomRectangleType&, size_t, 
                     FloatType&, FloatType&) const;

};


template<int D, typename KT, typename DT, typename QT, typename FT>
typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::eval(
    const GeomPointType &p, FloatType rtol) const {

  // initialize upper bound to be as if every point contributes K(0)
  // and lower bound to be as if every point contributes 0. 
  FloatType upper_bound = (data_tree_.root_)->attr_.weight(); 
  FloatType lower_bound = ConstantTraits<FloatType>::zero();

  // upper_bound and lower_bound are guaranteed to have the required tolerance
  single_tree(p, rtol, data_tree_.root_, data_tree_.bbox_, 0, upper_bound, lower_bound);

  // normalize the result properly
  FloatType result = lower_bound + (upper_bound - lower_bound) / 2;
  result *= KernelType::normalization;
  result /= (std::pow(bandwidth_, D) * data_tree_.size());

  return result;
}


// u and l are updated by amounts du and dl. it guarantees 
// exactly one of the following:
// 1. du and dl updates u and l such that the contribution of data points 
//    under r to the density is evaluated exactly. 
// 2. du and dl is such that |(u+du-l-dl)/(l+dl)| < rtol.
template<int D, typename KT, typename DT, typename QT, typename FT>
void KernelDensity<D,KT,DT,QT,FT>::single_tree(
    const GeomPointType &p, FloatType rtol, 
    DataNodeType *r, const GeomRectangleType &bounding_box, size_t depth,
    FloatType &u, FloatType &l) const {

  if (r == nullptr) { return; }

  GeomPointType proxy;

  // compute the upper and lower bound updates: these are 
  // the maximum and minimum possible contributions to the density
  // for any point in bounding_box
  proxy[0] = bounding_box.min_dist(p) / bandwidth_;
  FloatType max_val = kernel_.unnormalized_eval(proxy);
  FloatType du = r->attr_.weight()*(max_val-1);

  proxy[0] = bounding_box.max_dist(p) / bandwidth_;
  FloatType dl = r->attr_.weight()*kernel_.unnormalized_eval(proxy);

  // the updates are sufficient and there's no need to explore 
  // subtrees if
  // 1. the maximum contribution is less than epsilon. 
  // 2. the bound improvement already satisfies the required tolerance 
  if (max_val < std::numeric_limits<FloatType>::epsilon() ||
      std::abs((u+du-l-dl)/(l+dl)) < rtol) {
    u += du; l += dl; return;
  }

  // compute exact contributions if at a leaf, otherwise explore the
  // subtree recursively
  if (r->is_leaf()) {
    for (auto i = r->start_idx_; i <= r->end_idx_; ++i) {
      FloatType delta = kernel_.unnormalized_eval(
          (p - data_tree_.points_[i].point()) / bandwidth_
      );
      u+=delta; l+= delta;
    }
    u -= r->attr_.weight();
  } else {
    auto lower_bbox = bounding_box.lower_halfspace(depth, r->split_);
    auto upper_bbox = bounding_box.upper_halfspace(depth, r->split_);
    single_tree(p, rtol, r->left_, lower_bbox, (depth+1)%D, u, l);
    single_tree(p, rtol, r->right_, upper_bbox, (depth+1)%D, u, l);
  }

}

template<int D, typename KT, typename DT, typename QT, typename FT>
typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::naive_eval(const GeomPointType &p) const {
  FloatType total = ConstantTraits<FloatType>::zero();
  for (const auto &datum : data_tree_.points()) {
    total += kernel_.unnormalized_eval( (p - datum.point()) / bandwidth_  );
  }
  total *= KernelType::normalization;
  total /= (std::pow(bandwidth_, D) * data_tree_.size());
  return total;
}


template<int D, typename KT, typename DT, typename QT, typename FT>
inline typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::bandwidth() const { return bandwidth_; }

template<int D, typename KT, typename DT, typename QT, typename FT>
inline void KernelDensity<D,KT,DT,QT,FT>::set_bandwidth(FloatType bw) { 
  bandwidth_ = bw; 
}

template<int D, typename KT, typename DT, typename QT, typename FT>
void swap(KernelDensity<D,KT,DT,QT,FT> &lhs, KernelDensity<D,KT,DT,QT,FT> &rhs) {
  using std::swap;
  swap(lhs.bandwidth_, rhs.bandwidth_);
  swap(lhs.kernel_, rhs.kernel_);
  swap(lhs.data_tree_, rhs.data_tree_);
  return;
}

template<int D, typename KT, typename DT, typename QT, typename FT>
KernelDensity<D,KT,DT,QT,FT>::KernelDensity() : bandwidth_(1), kernel_(), data_tree_() {}

template<int D, typename KT, typename DT, typename QT, typename FT>
KernelDensity<D,KT,DT,QT,FT>::KernelDensity(
    const std::vector<DataPointType> &pts, FloatType bw, int leaf_max) 
  : bandwidth_(bw), kernel_(), data_tree_(pts, leaf_max) {}

template<int D, typename KT, typename DT, typename QT, typename FT>
KernelDensity<D,KT,DT,QT,FT>::KernelDensity(
    std::vector<DataPointType> &&pts, FloatType bw, int leaf_max) 
  : bandwidth_(bw), kernel_(), data_tree_(std::move(pts), leaf_max) {}

template<int D, typename KT, typename DT, typename QT, typename FT>
KernelDensity<D,KT,DT,QT,FT>::~KernelDensity() {}

}

#endif
