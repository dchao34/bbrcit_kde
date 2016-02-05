#ifndef BBRCITKDE_KERNELDENSITY_H__
#define BBRCITKDE_KERNELDENSITY_H__

#include <limits>
#include <queue>
#include <stack>
#include <tuple>
#include <utility>
#include <cassert>
#include <iostream>
#include <iomanip>

#include <Kdtree.h>
#include <DecoratedPoint.h>
#include <Attributes/KdeAttributes.h>
#include <EpanechnikovKernel.h>
#include <KdeTraits.h>
#include <FloatUtils.h>

namespace bbrcit {

template<int D, typename KT, typename AT, typename FT> class KernelDensity;

template<int D, typename KT, typename AT, typename FT>
void swap(KernelDensity<D,KT,AT,FT>&, KernelDensity<D,KT,AT,FT>&);

// KernelDensity<> implements a kernel density estimator in D dimensions 
// using a D dimensional Kdtree. 
template<int D, 
         typename KernelT=EpanechnikovKernel<D>,
         typename AttrT=KdeAttributes<double>,
         typename FloatT=double>
class KernelDensity {

  private:

    using KernelDensityType = KernelDensity<D,KernelT,AttrT,FloatT>;

    using KdtreeType = Kdtree<D,AttrT,FloatT>; 
    using TreeNodeType = typename KdtreeType::Node;

  public: 
    using DataPointType = typename KdtreeType::DataPointType;
    using GeomRectangleType = typename KdtreeType::RectangleType;
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

    // return the kde evaluated at point `p`. the result satisfies at least one of the following:
    // + the relative error is at most `rel_err`.
    // + the absolute error is at most `abs_err`.
    // otherwise, it will report to stderr that precision has been lost. 
    FloatT eval(const GeomPointType &p, FloatType rel_err, FloatType abs_err) const;

    // return the kde evaluate at points in `queries`. the relative error is at most `rel_err`.
    void eval(std::vector<DataPointType> &queries, FloatType rel_err) const;

    // return the kde naively evaluated at point `p`. slow... O(n^2); mainly for debugging
    FloatT naive_eval(const GeomPointType&) const;


  private:

    // bandwidth of the estimator
    FloatType bandwidth_;

    // the kernel function of the estimator 
    KernelType kernel_;

    // the kernel function of the estimator 
    KdtreeType data_tree_;

    // helper functions for the single point evaluation
    void single_tree(TreeNodeType*, const GeomPointType&, 
                     FloatType&, FloatType&, FloatType, FloatType, 
                     FloatType, FloatType) const;
    void single_tree_base(TreeNodeType*, const GeomPointType&,
                          FloatType, FloatType, FloatType&, FloatType&) const;
    void estimate_contributions(TreeNodeType*, const GeomPointType&, 
                                FloatType&, FloatType&) const;
    bool can_approximate(TreeNodeType*,
                         FloatType,FloatType,FloatType,FloatType,
                         FloatType,FloatType,FloatType,FloatType) const;
    void tighten_bounds(TreeNodeType*,
                        FloatType,FloatType,FloatType,FloatType,
                        FloatType&,FloatType&,FloatType,FloatType) const;

    // dual tree traversal for the multi point eval()
    void dual_tree(TreeNodeType*, TreeNodeType*, 
                    FloatType, FloatType, FloatType, 
                    const GeomRectangleType&, size_t, 
                    const GeomRectangleType&, size_t, 
                    KdtreeType&) const;

};

template<int D, typename KT, typename AT, typename FT>
inline typename KernelDensity<D,KT,AT,FT>::FloatType
KernelDensity<D,KT,AT,FT>::bandwidth() const { return bandwidth_; }

template<int D, typename KT, typename AT, typename FT>
inline void KernelDensity<D,KT,AT,FT>::set_bandwidth(FloatType bw) { 
  bandwidth_ = bw; 
}

template<int D, typename KT, typename AT, typename FT>
void swap(KernelDensity<D,KT,AT,FT> &lhs, KernelDensity<D,KT,AT,FT> &rhs) {
  using std::swap;
  swap(lhs.bandwidth_, rhs.bandwidth_);
  swap(lhs.kernel_, rhs.kernel_);
  swap(lhs.data_tree_, rhs.data_tree_);
  return;
}

template<int D, typename KT, typename AT, typename FT>
KernelDensity<D,KT,AT,FT>::KernelDensity() : bandwidth_(1), kernel_(), data_tree_() {}

template<int D, typename KT, typename AT, typename FT>
KernelDensity<D,KT,AT,FT>::KernelDensity(
    const std::vector<DataPointType> &pts, FloatType bw, int leaf_max) 
  : bandwidth_(bw), kernel_(), data_tree_(pts, leaf_max) {}

template<int D, typename KT, typename AT, typename FT>
KernelDensity<D,KT,AT,FT>::KernelDensity(
    std::vector<DataPointType> &&pts, FloatType bw, int leaf_max) 
  : bandwidth_(bw), kernel_(), data_tree_(std::move(pts), leaf_max) {}

template<int D, typename KT, typename AT, typename FT>
KernelDensity<D,KT,AT,FT>::~KernelDensity() {}


// single point kde evaluation. based on the algorithms given in the following:
// + ''Multiresolution Instance-Based Learning'' by Deng and Moore
// + ''Nonparametric Density Estimation: Toward Computational Tractability'' by Gray and Moore
template<int D, typename KT, typename AT, typename FT>
typename KernelDensity<D,KT,AT,FT>::FloatType
KernelDensity<D,KT,AT,FT>::eval(const GeomPointType &p, FloatType rel_err, FloatType abs_err) const {

  // initialization: 
  // + upper is such that all points contribute maximally (1)
  // + lower is such that all points contribute minimally (0)
  // (max is 1.0 because we do not incorporate the normalization until the end)
  FloatType upper = data_tree_.root_->attr_.weight();
  FloatType lower = ConstantTraits<FloatType>::zero();
  FloatType du = 1.0, dl = 0.0;

  // tighten the bounds by the single_tree algorithm
  single_tree(data_tree_.root_, p, upper, lower, du, dl, rel_err, abs_err);
  assert(lower <= upper); assert(lower >= 0);

  // take the mean of the bounds and remember to include the normalization
  FloatType normalization = KernelType::normalization / (std::pow(bandwidth_, D) * data_tree_.size());
  FloatType result = lower + (upper - lower) / 2;
  result *= normalization;

  // error reporting: notify the user of any loss of precision
  if (std::abs(normalization*(upper - lower)) > abs_err) {
    if (std::abs((upper-lower)/lower) > rel_err) {
      std::cerr << "Loss of relative precision with ";
      std::cerr << std::setprecision(5) << "relative error: " << std::abs((upper - lower)/lower);
      std::cerr << " (c.f. " << rel_err << ")";
      std::cerr << " when querying " << p << std::endl;
    } else {
      std::cerr << "Loss of absolute precision with ";
      std::cerr << std::setprecision(5) << "absolute error: " << std::abs(normalization*(upper - lower));
      std::cerr << " (c.f. " << abs_err << ")";
      std::cerr << " when querying " << p << std::endl;
    }
  }

  return result;

}

template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::single_tree(
    TreeNodeType *D_node, const GeomPointType &p, 
    FloatType &upper, FloatType &lower, 
    FloatType du, FloatType dl, 
    FloatType rel_err, FloatType abs_err) const {

  // update the individual contributions due to points in `D_node` towards the upper/lower 
  // bounds on the kde value at point `p`. 
  FloatType du_new, dl_new; 
  estimate_contributions(D_node, p, du_new, dl_new);

  // bound: approximate the total contribution due to `D_node` and decide whether to prune. 
  if (can_approximate(D_node, du_new, dl_new, du, dl, upper, lower, rel_err, abs_err)) { 

    // prune: still need to tighten the lower/upper bounds
    tighten_bounds(D_node, du_new, dl_new, du, dl, upper, lower, rel_err, abs_err);

    return; 
  }

  // branch: case 1: reached a leaf. brute force computation. 
  if (D_node->is_leaf()) {

    single_tree_base(D_node, p, du, dl, upper, lower);

  // branch: case 2: non-leaf. recursively tighten the bounds. 
  } else {

    // tighten the bounds for faster convergence
    tighten_bounds(D_node, du_new, dl_new, du, dl, upper, lower, rel_err, abs_err);

    // decide which halfspace is closer to the query
    TreeNodeType *closer = D_node->left_, *further = D_node->right_;
    if (D_node->left_->bbox_.min_dist(p) > D_node->right_->bbox_.min_dist(p)) {
      closer = D_node->right_; further = D_node->left_;
    }
    
    // recursively tighten the bounds, closer halfspace first 
    single_tree(closer, p, upper, lower, du_new, dl_new, rel_err, abs_err);
    single_tree(further, p, upper, lower, du_new, dl_new, rel_err, abs_err);

  }
}

template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::estimate_contributions(
    TreeNodeType *D_node, const GeomPointType &p, 
    FloatType &du, FloatType &dl) const {

  static GeomPointType proxy;

  // simply use the minimum(maximum) distance to the point to bound the 
  // maximum(minimum) contribution to the kernel evaluation
  proxy[0] = D_node->bbox_.min_dist(p) / bandwidth_;
  du = kernel_.unnormalized_eval(proxy);

  proxy[0] = D_node->bbox_.max_dist(p) / bandwidth_;
  dl = kernel_.unnormalized_eval(proxy);

}

// decide whether the current updates allow a prune
//
// + For the condition that gurantees the absolute errors, see 
//   Section 5. of ''Multiresolution Instance-Based Learning'' by Deng and Moore
//
// + For the condition that gurantees the relative errors, see 
//   Section 4.3. of ''Nonparametric Density Estimation: Toward Computational Tractability'' 
//   by Gray and Moore
template<int D, typename KT, typename AT, typename FT>
bool KernelDensity<D,KT,AT,FT>::can_approximate(
    TreeNodeType *D_node,
    FloatType du_new, FloatType dl_new, 
    FloatType du, FloatType dl, 
    FloatType upper, FloatType lower, 
    FloatType rel_err, FloatType abs_err) const {

  FloatType abs_tol = 2 * abs_err / data_tree_.size();

  // exclusion pruning guaranteeing that the absolute error <= abs_err
  if (std::abs(du_new) <= abs_tol) { return true; }

  // approximation pruning
  // condition 1: guarantee absolute error <= abs_err
  // condition 2: guarantee relative error <= rel_err
  if (std::abs(du_new - dl_new) <= abs_tol) { return true; }

  tighten_bounds(D_node, du_new, dl_new, du, dl, upper, lower, rel_err, abs_err);
  if (std::abs(upper - lower) <= abs_err || 
      std::abs(upper-lower) <= std::abs(lower)*rel_err) { return true; }

  return false;

}

// input invariants:
// + lower <= upper, dl <= du, dl_new <= du_new
// + dl <= dl_new, du >= du_new
//
// output invariants:
// + lower <= upper
template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::tighten_bounds(
    TreeNodeType *D_node,
    FloatType du_new, FloatType dl_new, 
    FloatType du, FloatType dl, 
    FloatType &upper, FloatType &lower, 
    FloatType rel_err, FloatType abs_err) const {

  // add the new contributions, but remember to subtract away the old ones
  lower += D_node->attr_.weight() * (dl_new - dl);
  upper += D_node->attr_.weight() * (du_new - du); 

  // the input invariants guarantee, mathematically, that lower <= upper.
  // however, roundoff error can cause approximate cancellations that break the gurantee.
  //
  // to enforce the output invariant, we explicitly set lower = upper 
  // if the cancellation overshoots. 
  if (lower > upper) { upper = lower; } 
}

// input invariants:
// + lower <= upper, dl <= du
//
// output invariants:
// + lower <= upper
template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::single_tree_base(
    TreeNodeType *D_node, const GeomPointType &p,
    FloatType du, FloatType dl, 
    FloatType &upper, FloatType &lower) const {

  FloatType delta;
  for (auto i = D_node->start_idx_; i <= D_node->end_idx_; ++i) {
    delta = kernel_.unnormalized_eval(
        (p - data_tree_.points_[i].point()) / bandwidth_
    );
    upper += delta; lower += delta;
  }
  upper -= D_node->attr_.weight() * du; lower -= D_node->attr_.weight() * dl;

  // see comment in tighten_bounds. 
  if (lower > upper) { upper = lower; }

}


// multi-point kde evaluation
template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::eval(
    std::vector<DataPointType> &queries, FloatType rel_err) const {

  for (auto &q : queries) { 
    auto attr = q.attributes();
    attr.set_lower(0); attr.set_upper(data_tree_.size());
    q.set_attributes(attr);
  }
  FloatType du = 1.0, dl = 0.0;

  KdtreeType query_tree(std::move(queries));

  dual_tree(query_tree.root_, data_tree_.root_, 
            du, dl, rel_err,
            query_tree.bbox_, 0, 
            data_tree_.bbox_, 0, 
            query_tree);

  FloatType normalization = KernelType::normalization; 
  normalization /= (std::pow(bandwidth_, D) * data_tree_.size());
  for (auto &q : query_tree.points_) { 
    auto attr = q.attributes();
    attr.set_lower(attr.lower()*normalization);
    attr.set_upper(attr.upper()*normalization);
    q.set_attributes(attr);
  }

  queries = std::move(query_tree.points_);

  return;

}


template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::dual_tree(
    TreeNodeType *Q_node, TreeNodeType *D_node, 
    FloatType du, FloatType dl, FloatType rel_tol, 
    const GeomRectangleType &Q_box, size_t Q_depth, 
    const GeomRectangleType &D_box, size_t D_depth, 
    KdtreeType &query_tree) const {

  static GeomPointType proxy;
  FloatType D_weight = D_node->attr_.weight();

  proxy[0] = D_box.min_dist(Q_box) / bandwidth_;
  FloatType du_new = kernel_.unnormalized_eval(proxy);

  proxy[0] = D_box.max_dist(Q_box) / bandwidth_;
  FloatType dl_new = kernel_.unnormalized_eval(proxy);

  Q_node->attr_.set_lower(Q_node->attr_.lower() + D_weight*(dl_new-dl));
  Q_node->attr_.set_upper(Q_node->attr_.upper() + D_weight*(du_new-du));

  // exclusion and approximation pruning
  if (almost_equal(du_new, ConstantTraits<FloatType>::zero()) ||
      std::abs(du_new - dl_new) < 2 * (Q_node->attr_.lower() + dl_new) * rel_tol / data_tree_.size()) {

    for (auto i = Q_node->start_idx_; i <= Q_node->end_idx_; ++i) {
      auto attr = query_tree.points_[i].attributes();
      attr.set_lower(attr.lower() + D_weight * dl_new);
      attr.set_upper(attr.upper() + D_weight * (du_new - 1));
      query_tree.points_[i].set_attributes(attr);
    }

    return;
  }

  // any node pair that reaches this point requires expansion 
  // to futher tighten their contributions

  // case 1: Q and D both leaves
  if (Q_node->is_leaf() && D_node->is_leaf()) {

    FloatType min_q = std::numeric_limits<FloatType>::max();
    FloatType max_q = std::numeric_limits<FloatType>::min();

    for (auto i = Q_node->start_idx_; i <= Q_node->end_idx_; ++i) {

      // update the contribution of each point in D
      auto attr = query_tree.points_[i].attributes();
      for (auto j = D_node->start_idx_; j <= D_node->end_idx_; ++j) {
        FloatType delta = kernel_.unnormalized_eval(
          (query_tree.points_[i].point() - data_tree_.points_[j].point()) / bandwidth_
        );
        attr.set_lower(attr.lower() + delta);
        attr.set_upper(attr.upper() + delta);
      }
      attr.set_upper(attr.upper() - D_weight);

      min_q = std::min(attr.lower(), min_q);
      max_q = std::max(attr.upper(), max_q);

      query_tree.points_[i].set_attributes(attr);

    }

    Q_node->attr_.set_lower(min_q);
    Q_node->attr_.set_upper(max_q);

  } else {

    // case 2: Q is leaf
    if (Q_node->is_leaf()) {

      GeomRectangleType D_lbox = D_box.lower_halfspace(D_depth, D_node->split_);
      GeomRectangleType D_ubox = D_box.upper_halfspace(D_depth, D_node->split_);

      // closer heuristic
      TreeNodeType *closer = D_node->left_; const GeomRectangleType *closer_r = &D_lbox;
      TreeNodeType *further = D_node->right_; const GeomRectangleType *further_r = &D_ubox;
      if (D_lbox.min_dist(Q_box) > D_ubox.min_dist(Q_box)) {
        closer = D_node->right_; closer_r = &D_ubox;
        further = D_node->left_; further_r = &D_lbox;
      }

      dual_tree(Q_node, closer, 
                du_new, dl_new, rel_tol, 
                Q_box, Q_depth, 
                *closer_r, (D_depth+1)%D, 
                query_tree);

      dual_tree(Q_node, further, 
                du_new, dl_new, rel_tol, 
                Q_box, Q_depth, 
                *further_r, (D_depth+1)%D, 
                query_tree);
    } else {

      Q_node->left_->attr_.set_lower(Q_node->attr_.lower());
      Q_node->left_->attr_.set_upper(Q_node->attr_.upper());
      Q_node->right_->attr_.set_lower(Q_node->attr_.lower());
      Q_node->right_->attr_.set_upper(Q_node->attr_.upper());

      GeomRectangleType Q_lbox = Q_box.lower_halfspace(Q_depth, Q_node->split_);
      GeomRectangleType Q_ubox = Q_box.upper_halfspace(Q_depth, Q_node->split_);

      // case 3: D is leaf
      if (D_node->is_leaf()) {

        dual_tree(Q_node->left_, D_node, 
                  du_new, dl_new, rel_tol, 
                  Q_lbox, (Q_depth+1)%D, 
                  D_box, D_depth, 
                  query_tree);

        dual_tree(Q_node->right_, D_node, 
                  du_new, dl_new, rel_tol, 
                  Q_ubox, (Q_depth+1)%D, 
                  D_box, D_depth, 
                  query_tree);

      // case 4: neither Q nor D are leaves
      } else {

        GeomRectangleType D_lbox = D_box.lower_halfspace(D_depth, D_node->split_);
        GeomRectangleType D_ubox = D_box.upper_halfspace(D_depth, D_node->split_);

        // Q left
        TreeNodeType *closer = D_node->left_; const GeomRectangleType *closer_r = &D_lbox;
        TreeNodeType *further = D_node->right_; const GeomRectangleType *further_r = &D_ubox;
        if (D_lbox.min_dist(Q_lbox) > D_ubox.min_dist(Q_lbox)) {
          closer = D_node->right_; closer_r = &D_ubox;
          further = D_node->left_; further_r = &D_lbox;
        }

        dual_tree(Q_node->left_, closer, 
                  du_new, dl_new, rel_tol, 
                  Q_lbox, (Q_depth+1)%D, 
                  *closer_r, (D_depth+1)%D, 
                  query_tree);

        dual_tree(Q_node->left_, further, 
                  du_new, dl_new, rel_tol, 
                  Q_lbox, (Q_depth+1)%D, 
                  *further_r, (D_depth+1)%D, 
                  query_tree);

        // Q right
        closer = D_node->left_; closer_r = &D_lbox;
        further = D_node->right_; further_r = &D_ubox;
        if (D_lbox.min_dist(Q_ubox) > D_ubox.min_dist(Q_ubox)) {
          closer = D_node->right_; closer_r = &D_ubox;
          further = D_node->left_; further_r = &D_lbox;
        }

        dual_tree(Q_node->right_, closer, 
                  du_new, dl_new, rel_tol, 
                  Q_ubox, (Q_depth+1)%D, 
                  *closer_r, (D_depth+1)%D, 
                  query_tree);

        dual_tree(Q_node->right_, further, 
                  du_new, dl_new, rel_tol, 
                  Q_ubox, (Q_depth+1)%D, 
                  *further_r, (D_depth+1)%D, 
                  query_tree);

      }

    }

  }

}


template<int D, typename KT, typename AT, typename FT>
typename KernelDensity<D,KT,AT,FT>::FloatType
KernelDensity<D,KT,AT,FT>::naive_eval(const GeomPointType &p) const {
  FloatType total = ConstantTraits<FloatType>::zero();
  for (const auto &datum : data_tree_.points()) {
    total += kernel_.unnormalized_eval( (p - datum.point()) / bandwidth_  );
  }
  total *= KernelType::normalization;
  total /= (std::pow(bandwidth_, D) * data_tree_.size());
  return total;
}


}

#endif
