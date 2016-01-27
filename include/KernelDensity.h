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
    FloatT eval2(const GeomPointType &p, FloatType rtol, size_t&) const;
    FloatT eval3(const GeomPointType &p, FloatType rtol, size_t&) const;
    FloatT eval4(const GeomPointType &p, FloatType rtol, size_t&) const;

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

    void walk() const;
    void walk_aux(DataNodeType*, const GeomRectangleType &, size_t, size_t&) const;

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

    void dual_tree(DataNodeType*, size_t, const GeomRectangleType&, 
                   QueryNodeType*, size_t, const GeomRectangleType&, 
                   QueryTreeType&, FloatType) const;

    struct QueueNode {
      DataNodeType *node_p;
      GeomRectangleType rectangle;
      size_t depth;
      double min_dist;
      double max_dist;
      bool operator<(const QueueNode &rhs) const {
        return rhs.min_dist < this->min_dist;
      }
    };

};

template<int D, typename KT, typename DT, typename QT, typename FT>
typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::eval4(
    const GeomPointType &p, FloatType rtol, size_t &cnt) const {

  cnt = 0;

  QueueNode q_node;
  q_node.node_p = data_tree_.root_;
  q_node.rectangle = data_tree_.bbox_;
  q_node.depth = 0;
  q_node.min_dist = data_tree_.bbox_.min_dist(p);
  q_node.max_dist = data_tree_.bbox_.max_dist(p);

  FloatType u = (data_tree_.root_)->attr_.weight();
  FloatType l = ConstantTraits<FloatType>::zero();

  std::priority_queue<QueueNode> q; q.push(q_node);
  while (!q.empty()) {

    DataNodeType *r = q.top().node_p; 
    GeomRectangleType bbox = q.top().rectangle; 
    size_t depth = q.top().depth; 
    FloatType min_d = q.top().min_dist; 
    FloatType max_d = q.top().max_dist; 
    q.pop();

    GeomPointType proxy;

    proxy[0] = min_d / bandwidth_;
    FloatType max_val = kernel_.unnormalized_eval(proxy);
    FloatType du_raw = r->attr_.weight() * max_val;
    FloatType du = du_raw - r->attr_.weight();

    proxy[0] = max_d / bandwidth_;
    FloatType dl = r->attr_.weight() * kernel_.unnormalized_eval(proxy);

    if (max_val < std::numeric_limits<FloatType>::epsilon() ||
        std::abs((du_raw-dl)/(l+dl)) < rtol) {
      u += du; l += dl; continue;
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
      auto lower_bbox = bbox.lower_halfspace(depth, r->split_);
      q_node.node_p = r->left_;
      q_node.rectangle = lower_bbox;
      q_node.depth = (depth+1)%D;
      q_node.min_dist = lower_bbox.min_dist(p);
      q_node.max_dist = lower_bbox.max_dist(p);
      q.push(q_node);
      auto upper_bbox = bbox.upper_halfspace(depth, r->split_);
      q_node.node_p = r->right_;
      q_node.rectangle = upper_bbox;
      q_node.depth = (depth+1)%D;
      q_node.min_dist = upper_bbox.min_dist(p);
      q_node.max_dist = upper_bbox.max_dist(p);
      q.push(q_node);
    }


  }

  std::cout << std::abs((u-l)/l) << std::endl;
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

    //single_tree3(p, rtol, r->left_, lower_bbox, (depth+1)%D, u, l, cnt);
    //single_tree3(p, rtol, r->right_, upper_bbox, (depth+1)%D, u, l, cnt);
  }

}

template<int D, typename KT, typename DT, typename QT, typename FT>
void KernelDensity<D,KT,DT,QT,FT>::walk() const {
  size_t cnt = 0;
  walk_aux(data_tree_.root_, data_tree_.bounding_box(), 0, cnt);
}

template<int D, typename KT, typename DT, typename QT, typename FT>
void KernelDensity<D,KT,DT,QT,FT>::walk_aux(
    DataNodeType *r, const GeomRectangleType &bounding_box, size_t depth, size_t &cnt) const {

  if (r->is_leaf()) { cnt += (r->end_idx_-r->start_idx_+1); return; }
  auto lower_bbox = bounding_box.lower_halfspace(depth, r->split_);
  auto upper_bbox = bounding_box.upper_halfspace(depth, r->split_);
  walk_aux(r->left_, lower_bbox, (depth+1)%D, cnt);
  walk_aux(r->right_, upper_bbox, (depth+1)%D, cnt);

}


template<int D, typename KT, typename DT, typename QT, typename FT>
typename KernelDensity<D,KT,DT,QT,FT>::FloatType
KernelDensity<D,KT,DT,QT,FT>::eval2(
    const GeomPointType &p, FloatType rtol, size_t &cnt) const {

  FloatType upper_bound = (data_tree_.root_)->attr_.weight(); 
  FloatType lower_bound = ConstantTraits<FloatType>::zero();

  cnt = 0;

  single_tree2(p, rtol, data_tree_.root_, data_tree_.bbox_, 0, upper_bound, lower_bound, 0, 1, cnt);

  FloatType result = lower_bound + (upper_bound - lower_bound) / 2;
  result *= KernelType::normalization;
  result /= (std::pow(bandwidth_, D) * data_tree_.size());

  return result;
}

template<int D, typename KT, typename DT, typename QT, typename FT>
void KernelDensity<D,KT,DT,QT,FT>::single_tree2(
    const GeomPointType &p, FloatType rtol, 
    DataNodeType *r, const GeomRectangleType &bounding_box, size_t depth,
    FloatType &u, FloatType &l, const FloatType &du, const FloatType &dl, size_t &cnt) const {

  if (r->is_leaf()) {
    for (auto i = r->start_idx_; i <= r->end_idx_; ++i) {
      FloatType delta = kernel_.unnormalized_eval(
          (p - data_tree_.points_[i].point()) / bandwidth_
      );
      u+=delta; l+= delta;
      ++cnt;
    }
    u -= r->attr_.weight();
    return;
  } 

  GeomPointType proxy;

  auto lower_bbox = bounding_box.lower_halfspace(depth, r->split_);
  proxy[0] = lower_bbox.min_dist(p) / bandwidth_;
  FloatType max_val_left = kernel_.unnormalized_eval(proxy);
  FloatType du_left = r->left_->attr_.weight()*(max_val_left-du);

  proxy[0] = lower_bbox.max_dist(p) / bandwidth_;
  FloatType min_val_left = kernel_.unnormalized_eval(proxy);
  FloatType dl_left = r->left_->attr_.weight()*(min_val_left-dl);

  u += du_left; l += dl_left;

  auto upper_bbox = bounding_box.upper_halfspace(depth, r->split_);
  proxy[0] = upper_bbox.min_dist(p) / bandwidth_;
  FloatType max_val_right = kernel_.unnormalized_eval(proxy);
  FloatType du_right = r->right_->attr_.weight()*(max_val_right-du);

  proxy[0] = upper_bbox.max_dist(p) / bandwidth_;
  FloatType min_val_right = kernel_.unnormalized_eval(proxy);
  FloatType dl_right = r->right_->attr_.weight()*(min_val_right-dl);

  u += du_right; l += dl_right;

  if (std::abs((u-l)/l) < rtol) { return; }

  if (max_val_left > std::numeric_limits<FloatType>::epsilon()) {
    single_tree2(p, rtol, r->left_, lower_bbox, (depth+1)%D, u, l, du_left, dl_left, cnt);
  }

  if (std::abs((u-l)/l) < rtol) { return; }

  if (max_val_right > std::numeric_limits<FloatType>::epsilon()) {
    single_tree2(p, rtol, r->right_, upper_bbox, (depth+1)%D, u, l, du_right, dl_right, cnt);
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


template<int D, typename KT, typename DT, typename QT, typename FT>
void KernelDensity<D,KT,DT,QT,FT>::eval(std::vector<QueryPointType> &query_points, FloatType rtol) const {

  FloatType upper_bound = (data_tree_.root_)->attr_.weight(); 
  FloatType lower_bound = ConstantTraits<FloatType>::zero();
  for (auto &p : query_points) { 
    p.set_attributes({p.attributes().weight(), lower_bound, upper_bound});
  }

  QueryTreeType query_tree(query_points);
  dual_tree(data_tree_.root_, 0, data_tree_.bounding_box(), 
            query_tree.root_, 0, query_tree.bounding_box(), 
            query_tree, rtol);

  FloatType normalization = KernelType::normalization / (std::pow(bandwidth_, D) * data_tree_.size());
  for (auto &p : query_points) { 
    p.set_attributes({p.attributes().weight(), 
                      p.attributes().lower() * normalization,
                      p.attributes().upper() * normalization});
  }

}

template<int D, typename KT, typename DT, typename QT, typename FT>
void KernelDensity<D,KT,DT,QT,FT>::dual_tree(
    DataNodeType *d_node, size_t d_depth, const GeomRectangleType &d_box,
    QueryNodeType *q_node, size_t q_depth, const GeomRectangleType &q_box,
    QueryTreeType &query_tree,
    FloatType rtol) const {

  GeomPointType proxy;
  proxy[0] = d_box.min_dist(q_box) / bandwidth_;
  FloatType max_val = kernel_.unnormalized_eval(proxy);
  FloatType du = d_node->attr_.weight()*(max_val-1);
  FloatType u = q_node->attr_.upper();

  proxy[0] = d_box.max_dist(q_box) / bandwidth_;
  FloatType dl = d_node->attr_.weight()*kernel_.unnormalized_eval(proxy);
  FloatType l = q_node->attr_.lower();

  if (max_val < std::numeric_limits<FloatType>::epsilon() ||
      std::abs((u+du-l-dl)/(l+dl)) < rtol) {
    for (auto i = q_node->start_idx_; i <= q_node->end_idx_; ++i) {
      auto &pt = query_tree.points_[i];
      pt.set_attributes({pt.attributes().weight(), 
                         pt.attributes().lower()+dl,
                         pt.attributes().upper()+du});
    }
    q_node->attr_.set_upper(u+du);
    q_node->attr_.set_lower(l+dl);
    return;
  }

  if (d_node->is_leaf() && q_node->is_leaf()) {
    for (auto i = q_node->start_idx_; i <= q_node->end_idx_; ++i) {
      FloatType delta = ConstantTraits<FloatType>::zero();
      auto &pt = query_tree.points_[i];
      for (auto j = d_node->start_idx_; j <= d_node->end_idx_; ++j) {
        delta += kernel_.unnormalized_eval(
            (pt.point() - data_tree_.points_[j].point()) / bandwidth_
        );
      }
      pt.set_attributes({
          pt.attributes().weight(), 
          pt.attributes().lower() + delta,
          pt.attributes().upper() + delta - d_node->attr_.weight()
      });
    }
  } else {

    if (!d_node->is_leaf() && !q_node->is_leaf()) {
      auto d_lbox = d_box.lower_halfspace(d_depth, d_node->split_);
      auto d_ubox = d_box.upper_halfspace(d_depth, d_node->split_);
      auto q_lbox = q_box.lower_halfspace(q_depth, q_node->split_);
      auto q_ubox = q_box.upper_halfspace(q_depth, q_node->split_);
      dual_tree(d_node->left_, (d_depth+1)%D, d_lbox, 
                q_node->left_, (q_depth+1)%D, q_lbox, 
                query_tree, rtol);
      dual_tree(d_node->left_, (d_depth+1)%D, d_lbox, 
                q_node->right_, (q_depth+1)%D, q_ubox, 
                query_tree, rtol);
      dual_tree(d_node->right_, (d_depth+1)%D, d_ubox, 
                q_node->left_, (q_depth+1)%D, q_lbox, 
                query_tree, rtol);
      dual_tree(d_node->right_, (d_depth+1)%D, d_ubox, 
                q_node->right_, (q_depth+1)%D, q_ubox, 
                query_tree, rtol);
    } else if (!q_node->is_leaf()) {
      auto q_lbox = q_box.lower_halfspace(q_depth, q_node->split_);
      auto q_ubox = q_box.upper_halfspace(q_depth, q_node->split_);
      dual_tree(d_node, d_depth, d_box, 
                q_node->left_, (q_depth+1)%D, q_lbox, 
                query_tree, rtol);
      dual_tree(d_node, d_depth, d_box, 
                q_node->right_, (q_depth+1)%D, q_ubox, 
                query_tree, rtol);
    } else {
      auto d_lbox = d_box.lower_halfspace(d_depth, d_node->split_);
      auto d_ubox = d_box.upper_halfspace(d_depth, d_node->split_);
      dual_tree(d_node->left_, (d_depth+1)%D, d_lbox, 
                q_node, q_depth, q_box, 
                query_tree, rtol);
      dual_tree(d_node->right_, (d_depth+1)%D, d_ubox, 
                q_node, q_depth, q_box, 
                query_tree, rtol);
    }

  }
}

}

#endif
