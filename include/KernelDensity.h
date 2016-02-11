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
#include <Kernels/EpanechnikovKernel.h>
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

    // construct a kernel density estimator with 
    // kernel `k` over the points `data`.
    KernelDensity(const std::vector<DataPointType> &data, int leaf_nmax=2);
    KernelDensity(std::vector<DataPointType> &&data, int leaf_nmax=2);

    // copy-control
    KernelDensity(const KernelDensityType&) = default;
    KernelDensity(KernelDensityType&&) noexcept = default;
    KernelDensityType& operator=(const KernelDensityType&) = default;
    KernelDensityType& operator=(KernelDensityType&&) = default;
    ~KernelDensity();

    // access/set the kernel
    const KernelType& kernel() const;
    KernelType& kernel();
    void set_kernel(const KernelType&);

    // returns the number of data points
    size_t size() const;

    // returns a const reference to the data points
    const std::vector<DataPointType>& points() const;

    // return the kde evaluated at point `p`. the result satisfies 
    // at least one of the following:
    //   + the relative error is at most `rel_err`.
    //   + the absolute error is at most `abs_err`.
    // otherwise, it will report to stderr that precision has been lost. 
    FloatT eval(DataPointType &p, 
                FloatType rel_err, FloatType abs_err) const;
    FloatT eval(const GeomPointType &p, 
                FloatType rel_err, FloatType abs_err) const;

    // return the kde evaluate at points in `queries`. 
    // the evaluataion at each point satisfies the guarantees 
    // as in the single point evaluation. 
    void eval(std::vector<DataPointType> &queries, 
              FloatType rel_err, FloatType abs_err) const;

    // return the kde naively evaluated at point `p`. 
    // slow... O(n^2); mainly for debugging
    FloatT naive_eval(DataPointType&) const;
    FloatT naive_eval(const GeomPointType&) const;


  private:

    // the kernel function of the estimator 
    KernelType kernel_;

    // the kernel function of the estimator 
    KdtreeType data_tree_;

    // helper functions for initialization
    void normalize_weights(std::vector<DataPointType>&);

    // helper functions for kde evaluataions
    void single_tree(const TreeNodeType*, const GeomPointType&, FloatType&, 
        FloatType&, FloatType, FloatType, FloatType, FloatType) const;
    void single_tree_base(const TreeNodeType*, const GeomPointType&,
        FloatType, FloatType, FloatType&, FloatType&) const;

    void dual_tree(const TreeNodeType*, TreeNodeType*, 
        FloatType, FloatType, FloatType, FloatType, KdtreeType&) const;
    void dual_tree_base(const TreeNodeType*, TreeNodeType*,
        FloatType, FloatType, KdtreeType&) const;

    bool can_approximate(const TreeNodeType*,
        FloatType,FloatType,FloatType,FloatType,
        FloatType,FloatType,FloatType,FloatType) const;
    bool can_approximate(const TreeNodeType*, const TreeNodeType*,
        FloatType,FloatType, FloatType, FloatType, 
        FloatType, FloatType) const;

    void tighten_bounds(const TreeNodeType*, FloatType, FloatType,
        FloatType,FloatType, FloatType&,FloatType&) const;
    void tighten_bounds(const TreeNodeType*, TreeNodeType*, 
        FloatType, FloatType, FloatType, FloatType) const;

    template<typename ObjT>
    void apply_closer_heuristic(
        TreeNodeType**, TreeNodeType**, const ObjT &) const;

    template<typename ObjT> 
    void estimate_contributions(const TreeNodeType*, const ObjT&, 
        FloatType&, FloatType&) const;

    void report_error(std::ostream&, const GeomPointType&,
        FloatType, FloatType, FloatType, FloatType) const;

};

template<int D, typename KT, typename AT, typename FT>
inline size_t KernelDensity<D,KT,AT,FT>::size() const { return data_tree_.size(); }

template<int D, typename KT, typename AT, typename FT>
inline const std::vector<typename KernelDensity<D,KT,AT,FT>::DataPointType>& 
KernelDensity<D,KT,AT,FT>::points() const {
  return data_tree_.points();
}

template<int D, typename KT, typename AT, typename FT>
void swap(KernelDensity<D,KT,AT,FT> &lhs, KernelDensity<D,KT,AT,FT> &rhs) {
  using std::swap;
  swap(lhs.kernel_, rhs.kernel_);
  swap(lhs.data_tree_, rhs.data_tree_);
  return;
}

template<int D, typename KT, typename AT, typename FT>
KernelDensity<D,KT,AT,FT>::KernelDensity() : 
  kernel_(), data_tree_() {}

template<int D, typename KT, typename AT, typename FT>
KernelDensity<D,KT,AT,FT>::KernelDensity(
    const std::vector<DataPointType> &pts, int leaf_max) 
  : kernel_() {

  std::vector<DataPointType> pts_normed = pts;
  normalize_weights(pts_normed);
  data_tree_ = KdtreeType(std::move(pts_normed), leaf_max);

}

template<int D, typename KT, typename AT, typename FT>
KernelDensity<D,KT,AT,FT>::KernelDensity(
    std::vector<DataPointType> &&pts, int leaf_max) 
  : kernel_() {

  normalize_weights(pts);
  data_tree_ = KdtreeType(std::move(pts), leaf_max);

}

template<int D, typename KT, typename AT, typename FT>
inline const typename KernelDensity<D,KT,AT,FT>::KernelType& 
KernelDensity<D,KT,AT,FT>::kernel() const {
  return kernel_;
}

template<int D, typename KT, typename AT, typename FT>
inline typename KernelDensity<D,KT,AT,FT>::KernelType& 
KernelDensity<D,KT,AT,FT>::kernel() {
  return const_cast<KernelType&>(
           static_cast<const KernelDensity<D,KT,AT,FT>&>(*this).kernel()
      );
}

template<int D, typename KT, typename AT, typename FT>
inline void KernelDensity<D,KT,AT,FT>::set_kernel(const KernelType &k) {
  kernel_ = k;
}

template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::normalize_weights(std::vector<DataPointType> &pts) {

  FloatType weight_total = ConstantTraits<FloatType>::zero();
  for (const auto &p : pts) { weight_total += p.attributes().weight(); }
  for (auto &p : pts) { 
    p.attributes().set_weight(
        p.attributes().weight() / weight_total
    );
  }

}

template<int D, typename KT, typename AT, typename FT>
KernelDensity<D,KT,AT,FT>::~KernelDensity() {}


template<int D, typename KT, typename AT, typename FT>
typename KernelDensity<D,KT,AT,FT>::FloatType
KernelDensity<D,KT,AT,FT>::eval(DataPointType &p, 
    FloatType rel_err, FloatType abs_err) const {

  FloatType result = eval(p.point(), rel_err, abs_err);
  p.attributes().set_upper(result); p.attributes().set_lower(result);
  return result;

}


// single point kde evaluation. based on the following algorithms:
// + ''Multiresolution Instance-Based Learning'' by Deng and Moore
// + ''Nonparametric Density Estimation: Toward Computational Tractability'' 
//   by Gray and Moore
template<int D, typename KT, typename AT, typename FT>
typename KernelDensity<D,KT,AT,FT>::FloatType
KernelDensity<D,KT,AT,FT>::eval(const GeomPointType &p, 
    FloatType rel_err, FloatType abs_err) const {

  // the contribution of some data point `d` to the kde at point `p`
  // is given by its weight scaled by the kernel value, which depends
  // on the distance between `d` and `p`. 
  //
  // (we will normalize the kernel once the computation is complete; 
  //  thus, we assume that the maximum value of any kernel is 1.0. )
  //
  // initialization: 
  // + upper: upper bound on the kde value. initially, every point 
  //          contributes its maximum weight. 
  // + lower: lower bound on the kde value. initially, every point 
  //          contributes none of its weight. 
  // + du: the upper bound on the proportion of weight each point contributes. 
  // + dl: the lower bound on the proportion of weight each point contributes. 
  FloatType upper = data_tree_.root_->attr_.weight();
  FloatType lower = ConstantTraits<FloatType>::zero();
  FloatType du = 1.0, dl = 0.0;

  // tighten the bounds by the single_tree algorithm. since we are computing
  // bounds before including the normalization, we need to scale abs_err accordingly
  FloatType normalization = kernel_.normalization();

  single_tree(data_tree_.root_, p, upper, lower, du, dl, 
              rel_err, abs_err / normalization);

  assert(lower <= upper); assert(lower >= 0);

  // take the mean of the bounds and remember to include the normalization
  FloatType result = normalization * (lower + (upper - lower) / 2);

  // error reporting: notify the user of any loss of precision
  report_error(std::cerr, p, 
               normalization*upper, 
               normalization*lower, 
               rel_err, abs_err);

  return result;

}

template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::single_tree(
    const TreeNodeType *D_node, const GeomPointType &p, 
    FloatType &upper, FloatType &lower, 
    FloatType du, FloatType dl, 
    FloatType rel_err, FloatType abs_err) const {

  // update the kernel contributions due to points in `D_node` 
  // towards the upper/lower bounds on the kde value at point `p`. 
  FloatType du_new, dl_new; 
  estimate_contributions(D_node, p, du_new, dl_new);

  // bound: approximate the total contribution due to `D_node` and 
  // decide whether to prune. 
  if (can_approximate(D_node, du_new, dl_new, du, dl, 
                      upper, lower, rel_err, abs_err)) { 

    // prune: still need to tighten the lower/upper bounds
    tighten_bounds(D_node, du_new, dl_new, du, dl, upper, lower);

    return; 
  }

  // branch: case 1: reached a leaf. brute force computation. 
  if (D_node->is_leaf()) {

    single_tree_base(D_node, p, du, dl, upper, lower);

  // branch: case 2: non-leaf. recursively tighten the bounds. 
  } else {

    // tighten the bounds for faster convergence
    tighten_bounds(D_node, du_new, dl_new, du, dl, upper, lower);

    // decide which halfspace is closer to the query
    TreeNodeType *closer = D_node->left_, *further = D_node->right_;
    apply_closer_heuristic(&closer, &further, p);
    
    // recursively tighten the bounds, closer halfspace first 
    single_tree(closer, p, upper, lower, du_new, dl_new, rel_err, abs_err);
    single_tree(further, p, upper, lower, du_new, dl_new, rel_err, abs_err);

  }
}

// input invariants:
// + lower <= upper, dl <= du
//
// output invariants:
// + lower <= upper
template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::single_tree_base(
    const TreeNodeType *D_node, const GeomPointType &p,
    FloatType du, FloatType dl, 
    FloatType &upper, FloatType &lower) const {

  FloatType delta;
  for (auto i = D_node->start_idx_; i <= D_node->end_idx_; ++i) {

    delta = kernel_.unnormalized_eval(p - data_tree_.points_[i].point());

    delta *= data_tree_.points_[i].attributes().weight();
    upper += delta; lower += delta;
  }
  upper -= D_node->attr_.weight() * du; 
  lower -= D_node->attr_.weight() * dl;

  // see comment in tighten_bounds. 
  if (lower > upper) { upper = lower; }

}


// multi-point kde evaluation
template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::eval(
    std::vector<DataPointType> &queries, 
    FloatType rel_err, FloatType abs_err) const {

  // initialize upper/lower bounds of individual queries to be
  // such that all data points contribute maximally/minimally
  for (auto &q : queries) { 
    q.attributes().set_lower(0);
    q.attributes().set_upper(data_tree_.root_->attr_.weight());
  }
  FloatType du = 1.0, dl = 0.0;

  // construct a query tree
  KdtreeType query_tree(std::move(queries));

  // tighten using the dual tree algorithm
  dual_tree(data_tree_.root_, query_tree.root_, 
            du, dl, rel_err, abs_err, query_tree);

  // remember to normalize before returning
  FloatType normalization = kernel_.normalization(); 
  for (auto &q : query_tree.points_) { 

    q.attributes().set_lower(q.attributes().lower()*normalization);
    q.attributes().set_upper(q.attributes().upper()*normalization);

    report_error(std::cerr, q.point(), 
                 q.attributes().upper(), q.attributes().lower(), 
                 rel_err, abs_err);
  }

  queries = std::move(query_tree.points_);

  return;
}


// tighten the contribution from all points in D_node to the upper/lower
// bounds of Q_node as well as each individual queries in Q_node
// `du`, `dl` are the present kernel contributions of D_node to bounds in Q_node. 
//
// the lower/upper bounds of Q_node is the min/max of all lower/upper 
// bounds of the individual queries 
template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::dual_tree(
    const TreeNodeType *D_node, TreeNodeType *Q_node, 
    FloatType du, FloatType dl, FloatType rel_err, FloatType abs_err,
    KdtreeType &query_tree) const {

  // update the kernel contributions due to D_node
  FloatType du_new, dl_new;
  estimate_contributions(D_node, Q_node->bbox_, du_new, dl_new);

  // BOUND: decide whether the approximation satsifies the error guarantees
  if (can_approximate(D_node, Q_node, 
        du_new, dl_new, du, dl, rel_err, abs_err)) {

    // tighten the lower/upper bound of Q_node itself
    tighten_bounds(D_node, Q_node, du_new, dl_new, du, dl);

    // tighten the individual queries
    double upper_q, lower_q;
    for (auto i = Q_node->start_idx_; i <= Q_node->end_idx_; ++i) {

      upper_q = query_tree.points_[i].attributes().upper();
      lower_q = query_tree.points_[i].attributes().lower();

      // du/dl are set to 1.0/0.0 because they were never 
      // updated since initialization
      tighten_bounds(D_node, du_new, dl_new, 1.0, 0.0, upper_q, lower_q);

      query_tree.points_[i].attributes().set_upper(upper_q);
      query_tree.points_[i].attributes().set_lower(lower_q);
    }

    return;
  }

  // BRANCH: any node pair that reaches this point requires expansion 
  // to further tighten their contributions

  // base case: Q and D both leaves; brute force
  if (Q_node->is_leaf() && D_node->is_leaf()) {

    dual_tree_base(D_node, Q_node, du, dl, query_tree);
    
  } else {

    // recursive cases below. 

    // case 1: Q is a leaf. tighten recursively with D_node's daughters.
    if (Q_node->is_leaf()) {

      // tighten Q_node bounds for faster convergence. 
      // this is just an optimization. 
      tighten_bounds(D_node, Q_node, du_new, dl_new, du, dl);

      // closer heuristic
      TreeNodeType *closer = D_node->left_, *further = D_node->right_;
      apply_closer_heuristic(&closer, &further, Q_node->bbox_);

      dual_tree(closer, Q_node, 
          du_new, dl_new, rel_err, abs_err, query_tree);
      dual_tree(further, Q_node, 
          du_new, dl_new, rel_err, abs_err, query_tree);

    } else {

      // in the cases below, proceed in two steps:
      //
      // + recursively tighten the contributions of D_node's daughters to 
      //   Q_node's daughters. 
      //
      // + obtain Q_node's bounds by taking the min/max daughter bounds. 

      // tighten bounds for faster convergence. this is just an optimization; 
      // one still needs to combine after recursion finishes.
      tighten_bounds(D_node, Q_node->left_, du_new, dl_new, du, dl);
      tighten_bounds(D_node, Q_node->right_, du_new, dl_new, du, dl);

      // case 2: D is a leaf
      if (D_node->is_leaf()) {

        dual_tree(D_node, Q_node->left_, 
            du_new, dl_new, rel_err, abs_err, query_tree);
        dual_tree(D_node, Q_node->right_, 
            du_new, dl_new, rel_err, abs_err, query_tree);

      // case 3: neither Q nor D are leaves
      } else {

        // tighten Q->left
        TreeNodeType *closer = D_node->left_, *further = D_node->right_;
        apply_closer_heuristic(&closer, &further, Q_node->left_->bbox_);

        dual_tree(closer, Q_node->left_, 
            du_new, dl_new, rel_err, abs_err, query_tree);
        dual_tree(further, Q_node->left_, 
            du_new, dl_new, rel_err, abs_err, query_tree);

        // tighten Q->right
        closer = D_node->left_; further = D_node->right_;
        apply_closer_heuristic(&closer, &further, Q_node->right_->bbox_);

        dual_tree(closer, Q_node->right_, 
            du_new, dl_new, rel_err, abs_err, query_tree);
        dual_tree(further, Q_node->right_, 
            du_new, dl_new, rel_err, abs_err, query_tree);

      }

      // combine the daughters' bounds to update Q_node's bounds
      Q_node->attr_.set_lower(
          std::min(Q_node->left_->attr_.lower(), 
                   Q_node->right_->attr_.lower()));
      Q_node->attr_.set_upper(
          std::max(Q_node->left_->attr_.upper(), 
                   Q_node->right_->attr_.upper()));
    }
  }
}


template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::dual_tree_base(
    const TreeNodeType *D_node, TreeNodeType *Q_node,
    FloatType du, FloatType dl, 
    KdtreeType &query_tree) const {

  FloatType min_q = std::numeric_limits<FloatType>::max();
  FloatType max_q = std::numeric_limits<FloatType>::min();

  FloatType lower_q, upper_q;
  for (auto i = Q_node->start_idx_; i <= Q_node->end_idx_; ++i) {

    // update the contribution of each point due to D_node
    upper_q = query_tree.points_[i].attributes().upper();
    lower_q = query_tree.points_[i].attributes().lower();

    single_tree_base(D_node, query_tree.points_[i].point(), 
        1.0, 0.0, upper_q, lower_q);

    query_tree.points_[i].attributes().set_lower(lower_q);
    query_tree.points_[i].attributes().set_upper(upper_q);

    min_q = std::min(lower_q, min_q);
    max_q = std::max(upper_q, max_q);

  }

  Q_node->attr_.set_lower(min_q);
  Q_node->attr_.set_upper(max_q);

}


template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::tighten_bounds(
    const TreeNodeType *D_node, TreeNodeType *Q_node,
    FloatType du_new, FloatType dl_new, 
    FloatType du, FloatType dl) const {

  FloatType upper = Q_node->attr_.upper();
  FloatType lower = Q_node->attr_.lower();

  tighten_bounds(D_node, du_new, dl_new, du, dl, upper, lower);

  Q_node->attr_.set_upper(upper);
  Q_node->attr_.set_lower(lower);
}


// input invariants:
// + lower <= upper, dl <= du, dl_new <= du_new
// + dl <= dl_new, du >= du_new
//
// output invariants:
// + lower <= upper
template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::tighten_bounds(
    const TreeNodeType *D_node,
    FloatType du_new, FloatType dl_new, 
    FloatType du, FloatType dl, 
    FloatType &upper, FloatType &lower) const {

  // add the new contributions, but remember to subtract away the old ones
  lower += D_node->attr_.weight() * (dl_new - dl);
  upper += D_node->attr_.weight() * (du_new - du); 

  // the input invariants guarantee, mathematically, that lower <= upper.
  // however, roundoff error (approx. cancellation) can break this gurantee.
  //
  // to enforce the output invariant, we explicitly set lower = upper 
  // if the cancellation overshoots. 
  if (lower > upper) { upper = lower; } 
}



template<int D, typename KT, typename AT, typename FT>
inline bool KernelDensity<D,KT,AT,FT>::can_approximate(
    const TreeNodeType *D_node, const TreeNodeType *Q_node,
    FloatType du_new, FloatType dl_new, 
    FloatType du, FloatType dl, 
    FloatType rel_err, FloatType abs_err) const {

  // safe to approximate only if all points can be approximated
  return can_approximate(D_node, du_new, dl_new, du, dl, 
                         Q_node->attr_.upper(), Q_node->attr_.lower(),
                         rel_err, abs_err);
}


// decide whether the current updates allow a prune
//
// + For the condition that gurantees the absolute errors, see 
//   Section 5. of Deng and Moore
//
// + For the condition that gurantees the relative errors, see 
//   Section 4.3. of Gray and Moore
template<int D, typename KT, typename AT, typename FT>
bool KernelDensity<D,KT,AT,FT>::can_approximate(
    const TreeNodeType *D_node,
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

  tighten_bounds(D_node, du_new, dl_new, du, dl, upper, lower);

  if (std::abs(upper-lower) <= abs_err || 
      std::abs(upper-lower) <= std::abs(lower)*rel_err) { return true; }

  return false;

}




template<int D, typename KT, typename AT, typename FT>
  template<typename ObjT>
inline void KernelDensity<D,KT,AT,FT>::apply_closer_heuristic(
    TreeNodeType **closer, TreeNodeType **further, const ObjT &obj) const {

  if ((*closer)->bbox_.min_dist(obj) > (*further)->bbox_.min_dist(obj)) {
    std::swap(*closer, *further);
  }

}


template<int D, typename KT, typename AT, typename FT>
  template<typename ObjT> 
void KernelDensity<D,KT,AT,FT>::estimate_contributions(
    const TreeNodeType *D_node, const ObjT &obj, 
    FloatType &du, FloatType &dl) const {

  static GeomPointType proxy;

  // use the minimum(maximum) distance to the argument in each 
  // dimension to bound the min/max kernel contributions

  for (int i = 0; i < D; ++i) { proxy[i] = D_node->bbox_.min_dist(i, obj); }
  du = kernel_.unnormalized_eval(proxy);

  for (int i = 0; i < D; ++i) { proxy[i] = D_node->bbox_.max_dist(i, obj); }
  dl = kernel_.unnormalized_eval(proxy);

}


template<int D, typename KT, typename AT, typename FT>
void KernelDensity<D,KT,AT,FT>::report_error(
    std::ostream &os, const GeomPointType &p,
    FloatType upper, FloatType lower, 
    FloatType rel_err, FloatType abs_err) const {

  if (std::abs(upper - lower) > abs_err) {
    if (lower) {
      if (std::abs((upper-lower)/lower) > rel_err) {
        os << std::setprecision(6);
        os << "Relative loss when querying " << p << ": " << std::endl;
        os << std::setprecision(15);
        os << "\tlower:   " << lower << std::endl;
        os << "\tupper:   " << upper << std::endl;
        os << "\tabs_err: " << std::abs(upper - lower) << " (c.f. " << abs_err << ")" << std::endl;
        os << "\trel_err: " << std::abs(upper - lower) / lower << " (c.f. " << rel_err << ")" << std::endl;
        os << std::endl;
      }
    } else {
      os << std::setprecision(6);
      os << "Absolute precision loss when querying " << p << ": " << std::endl;
      os << std::setprecision(15);
      os << "\tlower:   " << lower << std::endl;
      os << "\tupper:   " << upper << std::endl;
      os << "\tabs_err: " << std::abs(upper - lower) << " (c.f. " << abs_err << ")" << std::endl;
      os << "\trel_err: --- " << std::endl;
      os << std::endl;
    }
  }

}

template<int D, typename KT, typename AT, typename FT>
typename KernelDensity<D,KT,AT,FT>::FloatType
KernelDensity<D,KT,AT,FT>::naive_eval(DataPointType &p) const {
  FloatType result = naive_eval(p.point());
  p.attributes().set_upper(result);
  p.attributes().set_lower(result);
  return result;
}

template<int D, typename KT, typename AT, typename FT>
typename KernelDensity<D,KT,AT,FT>::FloatType
KernelDensity<D,KT,AT,FT>::naive_eval(const GeomPointType &p) const {
  FloatType total = ConstantTraits<FloatType>::zero();
  for (const auto &datum : data_tree_.points()) {
    total += 
      datum.attributes().weight() * 
      kernel_.unnormalized_eval( p - datum.point() );
  }
  total *= kernel_.normalization();
  return total;
}


}

#endif
