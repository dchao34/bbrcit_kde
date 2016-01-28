#ifndef BBRCITKDE_KERNELDENSITY_H__
#define BBRCITKDE_KERNELDENSITY_H__

#include <limits>
#include <queue>
#include <stack>
#include <tuple>
#include <utility>
#include <cassert>

#include <Kdtree.h>
#include <DecoratedPoint.h>
#include <PointWeights.h>
#include <QueryTreeAttributes.h>
#include <EpanechnikovKernel.h>
#include <KdeTraits.h>
#include <FloatUtils.h>

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

    using QueryTreeType = Kdtree<D,QAttrT,QAttrT,FloatT>; 
    using QueryNodeType = typename QueryTreeType::Node;

  public: 
    using DataPointType = typename DataTreeType::DataPointType;
    using QueryPointType = typename QueryTreeType::DataPointType;

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
    FloatT eval(const GeomPointType &p, FloatType rtol, size_t&) const;
    FloatT eval3(const GeomPointType &p, FloatType rtol, size_t&) const;
    FloatT eval4(const GeomPointType &p, FloatType rtol, size_t&) const;

    FloatT eval5(const GeomPointType &p, FloatType rtol, size_t&) const;
    FloatT eval6(const GeomPointType &p, FloatType rtol, size_t&) const;

    void eval(std::vector<QueryPointType>&, FloatType rtol) const;


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
                     FloatType&, FloatType&, size_t&) const;

    void single_tree2(const GeomPointType&, FloatType, 
                     DataNodeType*, const GeomRectangleType&, size_t, 
                     FloatType&, FloatType&, const FloatType&, const FloatType&, size_t&) const;

    void single_tree3(const GeomPointType&, FloatType, 
                     DataNodeType*, const GeomRectangleType&, size_t, 
                     FloatType&, FloatType&, size_t&) const;

    void single_tree_tighten(DataNodeType*, const GeomPointType&, 
                             FloatType&, FloatType&, FloatType, FloatType, 
                             FloatType, 
                             const GeomRectangleType&, size_t, size_t&) const;

    struct QueueNode {
      DataNodeType *p_dnode;
      GeomRectangleType rectangle;
      size_t depth;
      double min_dist;
      double max_dist;
      bool operator<(const QueueNode &rhs) const {
        return rhs.min_dist < this->min_dist;
      }
    };

    struct JobNode {
      DataNodeType *dnode_;
      double du_; 
      double dl_;
      GeomRectangleType bbox_;
      size_t depth_;
      double prio_;
      bool operator<(const JobNode &rhs) const {
        return rhs.prio_ < this->prio_;
      }
    };

};

template<int D, typename KT, typename DT, typename QT, typename FT>
typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::eval6(
    const GeomPointType &p, FloatType rtol, size_t &cnt) const {

  cnt = 0;

  FloatType upper = data_tree_.root_->attr_.weight();
  FloatType lower = ConstantTraits<FloatType>::zero();
  FloatType du = 1.0, dl = 0.0;

  std::priority_queue<JobNode> q; 
  FloatType max_prio = data_tree_.bbox_.max_dist(p);
  q.push({data_tree_.root_, du, dl, data_tree_.bbox_, 0, max_prio});

  GeomPointType proxy;
  while (!q.empty()) {

    JobNode curr = q.top(); q.pop();

    FloatType weight = curr.dnode_->attr_.weight();

    proxy[0] = curr.bbox_.min_dist(p) / bandwidth_;
    FloatType du_new = kernel_.unnormalized_eval(proxy);

    proxy[0] = curr.bbox_.max_dist(p) / bandwidth_;
    FloatType dl_new = kernel_.unnormalized_eval(proxy);

    upper += weight * (du_new - curr.du_);
    lower += weight * (dl_new - curr.dl_);

    // approximation pruning
    if (approxmately_equal(du_new, dl_new, rtol)) { continue; }

    // exclusion pruning
    if (almost_equal(du_new, ConstantTraits<FloatType>::zero())) { continue; }

    if (curr.dnode_->is_leaf()) {
      for (auto i = curr.dnode_->start_idx_; i <= curr.dnode_->end_idx_; ++i) {
        FloatType delta = kernel_.unnormalized_eval(
            (p - data_tree_.points_[i].point()) / bandwidth_
        );
        upper += delta; lower += delta;

        ++cnt;
      }
      upper -= weight * du_new; lower -= weight * dl_new;

    } else {

      GeomRectangleType lower_bbox = curr.bbox_.lower_halfspace(curr.depth_, curr.dnode_->split_);
      FloatType lower_prio = lower_bbox.min_dist(p);
      q.push({ curr.dnode_->left_, du_new, dl_new, lower_bbox, (curr.depth_+1)%D, lower_prio });

      GeomRectangleType upper_bbox = curr.bbox_.upper_halfspace(curr.depth_, curr.dnode_->split_);
      FloatType upper_prio = upper_bbox.min_dist(p);
      q.push({ curr.dnode_->right_, du_new, dl_new, upper_bbox, (curr.depth_+1)%D, upper_prio });
      
    }

  }

  FloatType result = lower + (upper - lower) / 2;
  result *= KernelType::normalization;
  result /= (std::pow(bandwidth_, D) * data_tree_.size());
  return result;

}


template<int D, typename KT, typename DT, typename QT, typename FT>
typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::eval5(
    const GeomPointType &p, FloatType rtol, size_t &cnt) const {

  cnt = 0;

  FloatType upper = data_tree_.root_->attr_.weight();
  FloatType lower = ConstantTraits<FloatType>::zero();
  FloatType du = 1.0, dl = 0.0;

  single_tree_tighten(data_tree_.root_, p, 
                      upper, lower, du, dl, rtol, 
                      data_tree_.bbox_, 0, 
                      cnt);

  FloatType result = lower + (upper - lower) / 2;
  result *= KernelType::normalization;
  result /= (std::pow(bandwidth_, D) * data_tree_.size());
  return result;

}

template<int D, typename KT, typename DT, typename QT, typename FT>
void KernelDensity<D,KT,DT,QT,FT>::single_tree_tighten(
    DataNodeType *r, const GeomPointType &p, 
    FloatType &upper, FloatType &lower, FloatType du, FloatType dl, 
    FloatType rtol, 
    const GeomRectangleType &bbox, size_t depth, 
    size_t &cnt) const {

  static GeomPointType proxy;
  FloatType weight = r->attr_.weight();

  proxy[0] = bbox.min_dist(p) / bandwidth_;
  FloatType du_new = kernel_.unnormalized_eval(proxy);

  proxy[0] = bbox.max_dist(p) / bandwidth_;
  FloatType dl_new = kernel_.unnormalized_eval(proxy);

  upper += weight * (du_new - du);
  lower += weight * (dl_new - dl);
  assert(du_new <= du);
  assert(dl_new >= dl);

  // approximation pruning
  if (approxmately_equal(du_new, dl_new, rtol)) { return; }

  // exclusion pruning
  if (almost_equal(du_new, ConstantTraits<FloatType>::zero())) { return; }

  if (r->is_leaf()) {
    for (auto i = r->start_idx_; i <= r->end_idx_; ++i) {
      FloatType delta = kernel_.unnormalized_eval(
          (p - data_tree_.points_[i].point()) / bandwidth_
      );
      upper += delta; lower += delta;
      ++cnt;
    }
    upper -= weight * du_new; lower -= weight * dl_new;
  } else {

    auto lower_bbox = bbox.lower_halfspace(depth, r->split_);
    auto upper_bbox = bbox.upper_halfspace(depth, r->split_);

    // closer heuristic
    DataNodeType *closer = r->left_; const GeomRectangleType *closer_r = &lower_bbox;
    DataNodeType *further = r->right_; const GeomRectangleType *further_r = &upper_bbox;
    if (lower_bbox.min_dist(p) > upper_bbox.min_dist(p)) {
      closer = r->right_; closer_r = &upper_bbox;
      further = r->left_; further_r = &lower_bbox;
    }
    single_tree_tighten(closer, p, 
                        upper, lower, du_new, dl_new, rtol, 
                        *closer_r, (depth+1)%D, cnt);
    single_tree_tighten(further, p, 
                        upper, lower, du_new, dl_new, rtol, 
                        *further_r, (depth+1)%D, cnt);

  }
}

template<int D, typename KT, typename DT, typename QT, typename FT>
typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::eval4(
    const GeomPointType &p, FloatType rtol, size_t &cnt) const {

  cnt = 0;

  FloatType u = (data_tree_.root_)->attr_.weight();
  FloatType l = ConstantTraits<FloatType>::zero();

  std::priority_queue<QueueNode> q; 
  q.push({data_tree_.root_, data_tree_.bbox_, 0, 
          data_tree_.bbox_.min_dist(p), 
          data_tree_.bbox_.max_dist(p)});

  while (!q.empty()) {

    QueueNode curr = q.top(); q.pop();
    FloatType weight = curr.p_dnode->attr_.weight();

    GeomPointType proxy;

    proxy[0] = curr.max_dist / bandwidth_;
    FloatType dl = weight * kernel_.unnormalized_eval(proxy);

    proxy[0] = curr.min_dist / bandwidth_;
    FloatType max_val = kernel_.unnormalized_eval(proxy);
    FloatType du_temp = weight*max_val;
    FloatType du = du_temp - weight;

    if (max_val < std::numeric_limits<FloatType>::epsilon() ||
        std::abs((du_temp-dl)/(l+dl)) < rtol) {
      u += du; l += dl; continue;
    }

    if (curr.p_dnode->is_leaf()) {
      for (auto i = curr.p_dnode->start_idx_; i <= curr.p_dnode->end_idx_; ++i) {
        FloatType delta = kernel_.unnormalized_eval(
            (p - data_tree_.points_[i].point()) / bandwidth_
        );
        u+=delta; l+= delta;
        ++cnt;
      }
      u -= weight;
    } else {

      auto lower_bbox = curr.rectangle.lower_halfspace(curr.depth, curr.p_dnode->split_);
      q.push({ curr.p_dnode->left_, lower_bbox, (curr.depth+1)%D, 
               lower_bbox.min_dist(p), lower_bbox.max_dist(p) });

      auto upper_bbox = curr.rectangle.upper_halfspace(curr.depth, curr.p_dnode->split_);
      q.push({ curr.p_dnode->right_, upper_bbox, (curr.depth+1)%D, 
               upper_bbox.min_dist(p), upper_bbox.max_dist(p) });
    }

  }

  FloatType result = l + (u - l) / 2;
  result *= KernelType::normalization;
  result /= (std::pow(bandwidth_, D) * data_tree_.size());
  return result;

}

template<int D, typename KT, typename DT, typename QT, typename FT>
typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::eval3(
    const GeomPointType &p, FloatType rtol, size_t &cnt) const {

  FloatType upper_bound = (data_tree_.root_)->attr_.weight(); 
  FloatType lower_bound = ConstantTraits<FloatType>::zero();

  cnt = 0;

  single_tree3(p, rtol, data_tree_.root_, data_tree_.bbox_, 0, upper_bound, lower_bound, cnt);

  FloatType result = lower_bound + (upper_bound - lower_bound) / 2;
  result *= KernelType::normalization;
  result /= (std::pow(bandwidth_, D) * data_tree_.size());

  return result;
}


// proof of concept for various pruning stages
template<int D, typename KT, typename DT, typename QT, typename FT>
void KernelDensity<D,KT,DT,QT,FT>::single_tree3(
    const GeomPointType &p, FloatType rtol, 
    DataNodeType *r, const GeomRectangleType &bounding_box, size_t depth,
    FloatType &u, FloatType &l, size_t &cnt) const {

  if (r == nullptr) { return; }

  GeomPointType proxy;

  proxy[0] = bounding_box.min_dist(p) / bandwidth_;
  FloatType max_val = kernel_.unnormalized_eval(proxy);
  FloatType du_temp = r->attr_.weight()*(max_val);
  FloatType du = r->attr_.weight()*(max_val-1);

  proxy[0] = bounding_box.max_dist(p) / bandwidth_;
  FloatType dl = r->attr_.weight()*kernel_.unnormalized_eval(proxy);

  // 1. constant mass pruning
  if (max_val < std::numeric_limits<FloatType>::epsilon() ||
      std::abs((du_temp-dl)/(l+dl)) < rtol) {
    u += du; l += dl; return;
  }

  if (r->is_leaf()) {
    for (auto i = r->start_idx_; i <= r->end_idx_; ++i) {
      FloatType delta = kernel_.unnormalized_eval(
          (p - data_tree_.points_[i].point()) / bandwidth_
      );
      u+=delta; l+= delta;
      ++cnt;
    }
    u -= r->attr_.weight();
  } else {
    auto lower_bbox = bounding_box.lower_halfspace(depth, r->split_);
    auto upper_bbox = bounding_box.upper_halfspace(depth, r->split_);

    // 2. closer heuristic
    DataNodeType *closer = r->left_; const GeomRectangleType *closer_r = &lower_bbox;
    DataNodeType *further = r->right_; const GeomRectangleType *further_r = &upper_bbox;
    if (lower_bbox.min_dist(p) > upper_bbox.min_dist(p)) {
      closer = r->right_; closer_r = &upper_bbox;
      further = r->left_; further_r = &lower_bbox;
    }
    single_tree3(p, rtol, closer, *closer_r, (depth+1)%D, u, l, cnt);
    single_tree3(p, rtol, further, *further_r, (depth+1)%D, u, l, cnt);

  }

}

template<int D, typename KT, typename DT, typename QT, typename FT>
typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::eval(
    const GeomPointType &p, FloatType rtol, size_t &cnt) const {

  // initialize upper bound to be as if every point contributes K(0)
  // and lower bound to be as if every point contributes 0. 
  FloatType upper_bound = (data_tree_.root_)->attr_.weight(); 
  FloatType lower_bound = ConstantTraits<FloatType>::zero();

  cnt = 0;

  // upper_bound and lower_bound are guaranteed to have the required tolerance
  single_tree(p, rtol, data_tree_.root_, data_tree_.bbox_, 0, upper_bound, lower_bound, cnt);

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
    FloatType &u, FloatType &l, size_t &cnt) const {

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
      ++cnt;
    }
    u -= r->attr_.weight();
  } else {
    auto lower_bbox = bounding_box.lower_halfspace(depth, r->split_);
    auto upper_bbox = bounding_box.upper_halfspace(depth, r->split_);
    single_tree(p, rtol, r->left_, lower_bbox, (depth+1)%D, u, l, cnt);
    single_tree(p, rtol, r->right_, upper_bbox, (depth+1)%D, u, l, cnt);
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
