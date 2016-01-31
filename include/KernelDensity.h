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
#include <KdeAttributes.h>
#include <EpanechnikovKernel.h>
#include <KdeTraits.h>
#include <FloatUtils.h>

namespace bbrcit {

template<int D, typename KT, typename DT, typename FT> class KernelDensity;

template<int D, typename KT, typename DT, typename FT>
void swap(KernelDensity<D,KT,DT,FT>&, KernelDensity<D,KT,DT,FT>&);

// KernelDensity<> implements a kernel density estimator in D dimensions using 
// a D dimensional Kdtree. 
template<int D, 
         typename KernelT=EpanechnikovKernel<D>,
         typename DAttrT=KdeAttributes<double>,
         typename FloatT=double>
class KernelDensity {

  private:

    using KernelDensityType = KernelDensity<D,KernelT,DAttrT,FloatT>;

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

    // evaluate the kde at point `p`. relative error is no more than `rel_err`,
    // and the relative cutoff tolerance is `rel_tol`.
    FloatT eval(const GeomPointType &p, FloatType rel_err, FloatType rel_tol, size_t&) const;


    // primarily for debugging:
    bool empty() const { return data_tree_.root_ == nullptr; }
    size_t size() const { return data_tree_.size(); }
    const std::vector<DataPointType>& data_points() const { return data_tree_.points(); }
    void report_leaves(std::vector<std::pair<size_t,size_t>> &l) const { data_tree_.report_leaves(l); }
    void range_search(
        const GeomRectangleType &q, 
        std::vector<DataPointType> &r) const { data_tree_.range_search(q, r); }
    const DAttrT& root_attributes() const { return data_tree_.root_->attr_; }

    // mainly for debugging: 
    // + naive kde evaluation; slow... O(n^2). 
    // + eval_recursive: recursive implementation of eval.
    FloatT naive_eval(const GeomPointType&) const;
    FloatT eval_recursive(const GeomPointType &p, FloatType rtol, size_t&) const;

  private:

    // bandwidth of the estimator
    FloatType bandwidth_;

    // the kernel function of the estimator 
    KernelType kernel_;

    // the kernel function of the estimator 
    DataTreeType data_tree_;


    // auxiliary function to eval_recursive()
    void single_tree_tighten(DataNodeType*, const GeomPointType&, 
                             FloatType&, FloatType&, FloatType, FloatType, 
                             FloatType, 
                             const GeomRectangleType&, size_t, size_t&) const;

    // priority queue node used in eval(). 
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

template<int D, typename KT, typename DT, typename FT>
typename KernelDensity<D,KT,DT,FT>::FloatType
KernelDensity<D,KT,DT,FT>::eval(const GeomPointType &p, 
                                   FloatType rel_err, FloatType rel_tol, size_t &cnt) const {

  cnt = 0;

  // upper and lower bounds on the value of f(p) (f is the kde)
  // du and dl are the contribution of every data point to the bounds. 
  FloatType upper = data_tree_.root_->attr_.weight();
  FloatType lower = ConstantTraits<FloatType>::zero();
  FloatType du = 1.0, dl = 0.0;

  // absolute error we are willing to make. the rationale follows the 'infernal zero' section of
  // https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
  FloatType abs_err = upper * rel_err;


  // priority queue for the data_tree nodes to processes. those that are closer 
  // to the query point p are higher prioriy. 
  // there are two types of nodes in the queue:
  // 1. points in the node still need their contributions to the bounds tightened. 
  // 2. points in the node that are done tightening according to rel_tol. 
  std::priority_queue<JobNode> q; 
  q.push({data_tree_.root_, du, dl, data_tree_.bbox_, 0, 0});

  // if the priority of a node is < max_prio, it is of type 1. otherwise it is type 2.
  FloatType max_prio = data_tree_.bbox_.max_dist(p) + 1;

  GeomPointType proxy;
  while (!q.empty()) {

    // if the the upper and lower bounds are 'close enough', meaning either 
    // in abs_err or in the relative error rel_err, then stop the loop. 
    if (approxmately_equal(lower, upper, rel_err, abs_err)) { break; }

    // process the highest prioiry node. two cases.
    JobNode curr = q.top(); q.pop();
    FloatType weight = curr.dnode_->attr_.weight();
    FloatType du_new, dl_new;


    // case 1: type 1 node. proceed to tighten its contribution to the bounds 
    if (curr.prio_ < max_prio) {

      // compute the new contributions
      proxy[0] = curr.bbox_.min_dist(p) / bandwidth_;
      du_new = kernel_.unnormalized_eval(proxy);

      proxy[0] = curr.bbox_.max_dist(p) / bandwidth_;
      dl_new = kernel_.unnormalized_eval(proxy);

      // update the bounds, but remember to subtract away its previous contribution
      upper += weight * (du_new - curr.du_);
      lower += weight * (dl_new - curr.dl_);

      // approximation and exclusion pruning: prune if the contributions is close enough. 
      // add max_prio to the current priority to indicate that it is now a type 2 node. 
      if (approxmately_equal(du_new, dl_new, rel_tol) || 
          almost_equal(du_new, ConstantTraits<FloatType>::zero())) { 
        q.push({ curr.dnode_, du_new, dl_new, curr.bbox_, curr.depth_, curr.prio_+max_prio });
        continue; 
      }

    // case 2: the contributions staty the same. proceed below. 
    } else {

      du_new = curr.du_; dl_new = curr.dl_;

    }

    // any node that reaches this point requires further tightening of their 
    // contributions to the upper and lower bounds. again, two cases. 

    // case 1: for a leaf, we have no choice but to do the naive computation. 
    // the contribution is now exact, and we do not insert it back into the queue. 
    if (curr.dnode_->is_leaf()) {
      for (auto i = curr.dnode_->start_idx_; i <= curr.dnode_->end_idx_; ++i) {
        FloatType delta = kernel_.unnormalized_eval(
            (p - data_tree_.points_[i].point()) / bandwidth_
        );
        upper += delta; lower += delta;

        ++cnt;
      }
      upper -= weight * du_new; lower -= weight * dl_new;

    // case 2: recursively tighten the children. 
    } else {

      GeomRectangleType lower_bbox = curr.bbox_.lower_halfspace(curr.depth_, curr.dnode_->split_);
      FloatType lower_prio = lower_bbox.min_dist(p);
      q.push({ curr.dnode_->left_, du_new, dl_new, lower_bbox, (curr.depth_+1)%D, lower_prio });

      GeomRectangleType upper_bbox = curr.bbox_.upper_halfspace(curr.depth_, curr.dnode_->split_);
      FloatType upper_prio = upper_bbox.min_dist(p);
      q.push({ curr.dnode_->right_, du_new, dl_new, upper_bbox, (curr.depth_+1)%D, upper_prio });
      
    }

  }

  // return the mean of the upper and lower bounds 
  FloatType result = lower + (upper - lower) / 2;
  result *= KernelType::normalization;
  result /= (std::pow(bandwidth_, D) * data_tree_.size());

  return result;

}


// recursive implementation of eval(). mainly for cross checks. 
template<int D, typename KT, typename DT, typename FT>
typename KernelDensity<D,KT,DT,FT>::FloatType
KernelDensity<D,KT,DT,FT>::eval_recursive(
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

template<int D, typename KT, typename DT, typename FT>
void KernelDensity<D,KT,DT,FT>::single_tree_tighten(
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

template<int D, typename KT, typename DT, typename FT>
typename KernelDensity<D,KT,DT,FT>::FloatType
KernelDensity<D,KT,DT,FT>::naive_eval(const GeomPointType &p) const {
  FloatType total = ConstantTraits<FloatType>::zero();
  for (const auto &datum : data_tree_.points()) {
    total += kernel_.unnormalized_eval( (p - datum.point()) / bandwidth_  );
  }
  total *= KernelType::normalization;
  total /= (std::pow(bandwidth_, D) * data_tree_.size());
  return total;
}


template<int D, typename KT, typename DT, typename FT>
inline typename KernelDensity<D,KT,DT,FT>::FloatType
KernelDensity<D,KT,DT,FT>::bandwidth() const { return bandwidth_; }

template<int D, typename KT, typename DT, typename FT>
inline void KernelDensity<D,KT,DT,FT>::set_bandwidth(FloatType bw) { 
  bandwidth_ = bw; 
}

template<int D, typename KT, typename DT, typename FT>
void swap(KernelDensity<D,KT,DT,FT> &lhs, KernelDensity<D,KT,DT,FT> &rhs) {
  using std::swap;
  swap(lhs.bandwidth_, rhs.bandwidth_);
  swap(lhs.kernel_, rhs.kernel_);
  swap(lhs.data_tree_, rhs.data_tree_);
  return;
}

template<int D, typename KT, typename DT, typename FT>
KernelDensity<D,KT,DT,FT>::KernelDensity() : bandwidth_(1), kernel_(), data_tree_() {}

template<int D, typename KT, typename DT, typename FT>
KernelDensity<D,KT,DT,FT>::KernelDensity(
    const std::vector<DataPointType> &pts, FloatType bw, int leaf_max) 
  : bandwidth_(bw), kernel_(), data_tree_(pts, leaf_max) {}

template<int D, typename KT, typename DT, typename FT>
KernelDensity<D,KT,DT,FT>::KernelDensity(
    std::vector<DataPointType> &&pts, FloatType bw, int leaf_max) 
  : bandwidth_(bw), kernel_(), data_tree_(std::move(pts), leaf_max) {}

template<int D, typename KT, typename DT, typename FT>
KernelDensity<D,KT,DT,FT>::~KernelDensity() {}


}

#endif
