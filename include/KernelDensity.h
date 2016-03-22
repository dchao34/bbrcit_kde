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
#include <Attributes/AdaKdeAttributes.h>
#include <Kernels/EpanechnikovKernel.h>
#include <Kernels/KernelTraits.h>
#include <Kernels/ConvKernelAssociator.h>
#include <KdeTraits.h>
#include <FloatUtils.h>

#ifdef __CUDACC__
#include <CudaDirectKde.h>
#endif

namespace bbrcit {

template<int D, typename KT, typename FT, typename AT> class KernelDensity;

template<int D, typename KT, typename FT, typename AT>
void swap(KernelDensity<D,KT,FT,AT>&, KernelDensity<D,KT,FT,AT>&);

// KernelDensity<> implements a kernel density estimator in D dimensions 
// using a D dimensional Kdtree. 
template<int D, 
         typename KernelT=EpanechnikovKernel<D,double>,
         typename FloatT=double,
         typename AttrT=AdaKdeAttributes<FloatT>>
class KernelDensity {

  public: 

    using FloatType = FloatT;
    using KdtreeType = Kdtree<D,AttrT,FloatT>; 
    using DataPointType = typename KdtreeType::DataPointType;
    using KernelType = KernelT;
    using KernelFloatType = typename KernelType::FloatType;

    static constexpr int dim() { return D; }

  private:

    using KernelDensityType = KernelDensity<D,KernelT,FloatT,AttrT>;
    using TreeNodeType = typename KdtreeType::Node;
    using GeomPointType = typename KdtreeType::DataPointType::PointType;

    friend void swap<>(KernelDensityType&, KernelDensityType&);

  public:

    // default constructor: estimator on an empty set. 
    KernelDensity();

    // construct a kernel density estimator with 
    // kernel `k` over the points `data`.
    KernelDensity(std::vector<DataPointType> data, int leaf_nmax=2);
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


    // return the direct kde evaluated at point `p` or at points in `queries`
    FloatT direct_eval(DataPointType &p) const;
#ifndef __CUDACC__
    void direct_eval(std::vector<DataPointType> &queries) const;
#else
    void direct_eval(std::vector<DataPointType> &queries, size_t block_size=128) const;
#endif


    // return the kde evaluated at point `p` or at points in `queries`.
    // the result satisfies at least one of the following:
    //   + the relative error is at most `rel_err`.
    //   + the absolute error is at most `abs_err`.
    // otherwise, it will report to stderr that precision has been lost. 
    //
    // Note: for multi point queries, the method taking a vector<> first constructs 
    // KdtreeType<> before evaluation. Since such a construction uses randomized 
    // partitioning, the exact structure of the tree varies between calls. For 
    // uses where the exact structure matters, use the overloaded method taking a
    // KdtreeType<> argrument. 
    FloatT eval(DataPointType &p, FloatType rel_err, FloatType abs_err) const;
#ifndef __CUDACC__
    void eval(std::vector<DataPointType> &queries, 
              FloatType rel_err, FloatType abs_err, int leaf_nmax=2) const;
    void eval(KdtreeType &query_tree, FloatType rel_err, FloatType abs_err) const;
#else
    void eval(std::vector<DataPointType> &queries, 
              FloatType rel_err, FloatType abs_err, int leaf_nmax=2, 
              size_t block_size=128) const;
    void eval(KdtreeType &query_tree, FloatType rel_err, FloatType abs_err, 
              size_t block_size=128) const;
#endif


    // convert this object into an adaptive kernel. this method should be called 
    // once and is not reversible. 
#ifndef __CUDACC__ 
    void adapt_density(FloatType alpha, 
        FloatType rel_err=1e-6, FloatType abs_err=1e-6);
#else
    void adapt_density(FloatType alpha, 
        FloatType rel_err=1e-6, FloatType abs_err=1e-6, size_t block_size=128);
#endif


    // compute the cross validation score for the current 
    // kernel configuration. likelihood and least squares are available. 
#ifndef __CUDACC__ 
    FloatType likelihood_cross_validate(
        FloatType rel_err=1e-6, FloatType abs_err=1e-6) const;
    FloatType leastsquares_cross_validate(
        FloatType rel_err=1e-6, FloatType abs_err=1e-6) const;
#else
    FloatType likelihood_cross_validate(
        FloatType rel_err=1e-6, FloatType abs_err=1e-6, size_t block_size=128) const;
    FloatType leastsquares_cross_validate(
        FloatType rel_err=1e-6, FloatType abs_err=1e-6, size_t block_size=128) const;
#endif


  private:

    // the kernel function of the estimator 
    KernelType kernel_;

    // the kernel function of the estimator 
    KdtreeType data_tree_;

    // helper functions for initialization
    void initialize_attributes(std::vector<DataPointType>&);
    void normalize_weights(std::vector<DataPointType>&);


    // helper functions for direct kde evaluations
    // ------------------------------------------

    template<typename KernT> 
      FloatT direct_eval(const GeomPointType&, const KernT&) const;

#ifndef __CUDACC__
    template<typename KernT>
      void direct_eval(std::vector<DataPointType>&, const KernT&) const;
#else
    template<typename KernT>
      void direct_eval(std::vector<DataPointType>&, const KernT&, size_t) const;
#endif

    // helper functions for tree kde evaluataions
    // ------------------------------------------

    // single tree

    template <typename KernT>
      FloatT eval(const GeomPointType&, const KernT&, FloatType, FloatType) const;

    template <typename KernT>
      void single_tree(
          const TreeNodeType*, const GeomPointType&, const KernT&, 
          FloatType&, FloatType&, FloatType, FloatType, FloatType, FloatType) const;

    template <typename KernT>
      void single_tree_base(
          const TreeNodeType*, const GeomPointType&, const KernT&,
          FloatType, FloatType, FloatType&, FloatType&) const;


    // dual tree
#ifndef __CUDACC__ 
    template<typename KernT>
      void eval(KdtreeType&, const KernT&, FloatType, FloatType) const;

    template<typename KernT>
      void dual_tree(const TreeNodeType*, TreeNodeType*, const KernT&,
          FloatType, FloatType, FloatType, FloatType, KdtreeType&) const;

    template<typename KernT>
      void dual_tree_base(const TreeNodeType*, TreeNodeType*, const KernT&,
          FloatType, FloatType, KdtreeType&) const;
#else

    template<typename KernT>
      void eval(KdtreeType&, const KernT&, FloatType, FloatType, size_t) const;

    template<typename KernT>
      void dual_tree(const TreeNodeType*, TreeNodeType*, const KernT&,
          FloatType, FloatType, FloatType, FloatType, KdtreeType&, 
          CudaDirectKde<D,KernelFloatType,KernT>&, 
          std::vector<KernelFloatType>&,size_t) const;

    template<typename KernT>
      void dual_tree_base(const TreeNodeType*, TreeNodeType*, const KernT&,
          FloatType, FloatType, KdtreeType&, 
          CudaDirectKde<D,KernelFloatType,KernT>&, 
          std::vector<KernelFloatType>&,size_t) const;
#endif

    // general
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

    template<typename ObjT, typename KernT> 
    void estimate_contributions(
        const TreeNodeType*, const ObjT&, const KernT&,
        FloatType&, FloatType&) const;

    void report_error(std::ostream&, const GeomPointType&,
        FloatType, FloatType, FloatType, FloatType) const;


};

// perform least squares cross validation on the current kernel
// configuration. 
// NOTE: the weights of the data points are assumed to be normalized. 
// In particular, this holds before adapt_density() is called. Perhaps 
// consider removing this caveat? 
template<int D, typename KT, typename FT, typename AT>
typename KernelDensity<D,KT,FT,AT>::FloatType
KernelDensity<D,KT,FT,AT>::leastsquares_cross_validate( 
#ifndef __CUDACC__ 
    FloatType rel_err, FloatType abs_err
#else
    FloatType rel_err, FloatType abs_err, size_t block_size
#endif
    ) const {


  // construct a query tree out of the data tree to perform 
  // dual tree self-evaluation. we use copy-assignment since we would like
  // to preserve the same ordering of points in both trees. 
  KdtreeType query_tree = data_tree_;

  // compute the leave one out contribution
  // --------------------------------------

  // all pairs computation using the default kernel
#ifndef __CUDACC__
  eval(query_tree, kernel_, rel_err, abs_err);
#else
  eval(query_tree, kernel_, rel_err, abs_err, block_size);
#endif

  FloatType val = 0.0;

  // compute leave one out score
  FloatType llo_cv = ConstantTraits<FloatType>::zero();
  for (size_t i = 0; i < query_tree.points_.size(); ++i) {

    // the dual tree gives contributions from all points; must 
    // subtract away the self contribution
    val = query_tree.points_[i].attributes().middle();
    val -= data_tree_.points_[i].attributes().weight() * kernel_.normalization();

    // contribution is weighted
    llo_cv += data_tree_.points_[i].attributes().weight() * val;
  }

  // compute the square integral contribution
  // ----------------------------------------

  // induce the convolution kernel out of the default kernel
  typename ConvKernelAssociator<KernelType>::ConvKernelType conv_kernel = 
    ConvKernelAssociator<KernelType>::make_convolution_kernel(kernel_);

  // all pairs computation using the convolution kernel
#ifndef __CUDACC__
  eval(query_tree, conv_kernel, rel_err, abs_err);
#else
  eval(query_tree, conv_kernel, rel_err, abs_err, block_size);
#endif

  // compute square integral score
  FloatType sq_cv = ConstantTraits<FloatType>::zero();
  for (size_t i = 0; i < query_tree.points_.size(); ++i) {

    val = query_tree.points_[i].attributes().middle();

    // contribution is weighted
    sq_cv += data_tree_.points_[i].attributes().weight() * val;
  }

  return sq_cv - 2*llo_cv;

}



// perform likelihood cross validation on the current kernel
// configuration. 
// NOTE: the weights of the data points are assumed to be normalized. 
// In particular, this holds before adapt_density() is called. Perhaps 
// consider removing this caveat? 
template<int D, typename KT, typename FT, typename AT>
typename KernelDensity<D,KT,FT,AT>::FloatType
KernelDensity<D,KT,FT,AT>::likelihood_cross_validate( 
#ifndef __CUDACC__ 
    FloatType rel_err, FloatType abs_err
#else
    FloatType rel_err, FloatType abs_err, size_t block_size
#endif
    ) const {

  // construct a query tree out of the data tree to perform 
  // dual tree self-evaluation. we use copy-assignment since we would like
  // to preserve the same ordering of points in both trees. 
  KdtreeType query_tree = data_tree_;

#ifndef __CUDACC__
  eval(query_tree, kernel_, rel_err, abs_err);
#else
  eval(query_tree, kernel_, rel_err, abs_err, block_size);
#endif

  // compute the cross validation score
  FloatType cv = ConstantTraits<FloatType>::zero();

  FloatType cv_i;
  for (size_t i = 0; i < query_tree.points_.size(); ++i) {

    // the dual tree gives contributions from all points; must 
    // subtract away the self contribution
    cv_i = query_tree.points_[i].attributes().middle();
    cv_i -= data_tree_.points_[i].attributes().weight() * kernel_.normalization();

    // the cross validation score is the log of the leave one out contribution
    cv += data_tree_.points_[i].attributes().weight() * std::log(cv_i);
  }

  return cv;
}



// Calling this method repurposes `this` KernelDensity object to become 
// an adaptive kernel density. In particular, the following attributes in
// the data_tree must be updated:
//
// + For each node, update the min/max local bandwidth corrections of 
//   points under it. 
//
// + Update weights for each point and each node. Node weights are 
//   induced from point weights, while point weights are adjusted by 
//   scaling. e.g. if the `i`th point currently has weight `w_i` and 
//   local bandwidth correction `abw_i`, then update the weight to 
//   `w_i / pow(abw_i,D)`. 
//
// This prescription is described in page 101 of Silverman's book
// `Density Estimation for Statistics and Data Analysis`. 
template<int D, typename KT, typename FT, typename AT>
void KernelDensity<D,KT,FT,AT>::adapt_density(
#ifndef __CUDACC__
    FloatType alpha, FloatType rel_err, FloatType abs_err
#else
    FloatType alpha, FloatType rel_err, FloatType abs_err,
    size_t block_size
#endif
    ) {

  // compute pilot estimate
  // ----------------------

  // construct a query tree out of the data tree to perform 
  // dual tree self-evaluation. we use copy-assignment since we would like
  // to preserve the same ordering of points in both trees. 
  KdtreeType query_tree = data_tree_;

#ifndef __CUDACC__
  eval(query_tree, kernel_, rel_err, abs_err);
#else
  eval(query_tree, kernel_, rel_err, abs_err, block_size);
#endif

  // compute local bandwidth corrections
  // -----------------------------------

  FloatType g = 0;
  std::vector<FloatType> local_bw(query_tree.points_.size());
  for (size_t i = 0; i < query_tree.points_.size(); ++i) {
    local_bw[i] = query_tree.points_[i].attributes().middle();
    g += data_tree_.points_[i].attributes().weight() * std::log(local_bw[i]);
  }
  g = std::exp(g);

  for (auto &bw : local_bw) {
    bw = std::pow(bw/g, -alpha);
  }

  // update data tree attributes
  // ---------------------------
  for (size_t i = 0; i < data_tree_.points_.size(); ++i) {

    // local bandwidth corrections
    data_tree_.points_[i].attributes().set_lower_abw(local_bw[i]);
    data_tree_.points_[i].attributes().set_upper_abw(local_bw[i]);

    // scale wieghts
    data_tree_.points_[i].attributes().set_weight(
      data_tree_.points_[i].attributes().weight() * pow(local_bw[i], -D));
  }

  // update node attributes
  data_tree_.refresh_node_attributes(data_tree_.root_);


  return;
}

template<int D, typename KT, typename FT, typename AT>
inline size_t KernelDensity<D,KT,FT,AT>::size() const { return data_tree_.size(); }

template<int D, typename KT, typename FT, typename AT>
inline const std::vector<typename KernelDensity<D,KT,FT,AT>::DataPointType>& 
KernelDensity<D,KT,FT,AT>::points() const {
  return data_tree_.points();
}

template<int D, typename KT, typename FT, typename AT>
void swap(KernelDensity<D,KT,FT,AT> &lhs, KernelDensity<D,KT,FT,AT> &rhs) {
  using std::swap;
  swap(lhs.kernel_, rhs.kernel_);
  swap(lhs.data_tree_, rhs.data_tree_);
  return;
}

template<int D, typename KT, typename FT, typename AT>
KernelDensity<D,KT,FT,AT>::KernelDensity() : 
  kernel_(), data_tree_() {}

template<int D, typename KT, typename FT, typename AT>
KernelDensity<D,KT,FT,AT>::KernelDensity(
    std::vector<DataPointType> pts, int leaf_max) 
  : kernel_() {

  initialize_attributes(pts);
  data_tree_ = KdtreeType(std::move(pts), leaf_max);

}

template<int D, typename KT, typename FT, typename AT>
KernelDensity<D,KT,FT,AT>::KernelDensity(
    std::vector<DataPointType> &&pts, int leaf_max) 
  : kernel_() {

  initialize_attributes(pts);
  data_tree_ = KdtreeType(std::move(pts), leaf_max);

}

template<int D, typename KT, typename FT, typename AT>
inline const typename KernelDensity<D,KT,FT,AT>::KernelType& 
KernelDensity<D,KT,FT,AT>::kernel() const {
  return kernel_;
}

template<int D, typename KT, typename FT, typename AT>
inline typename KernelDensity<D,KT,FT,AT>::KernelType& 
KernelDensity<D,KT,FT,AT>::kernel() {
  return const_cast<KernelType&>(
           static_cast<const KernelDensity<D,KT,FT,AT>&>(*this).kernel()
      );
}

template<int D, typename KT, typename FT, typename AT>
inline void KernelDensity<D,KT,FT,AT>::set_kernel(const KernelType &k) {
  kernel_ = k;
}

template<int D, typename KT, typename FT, typename AT>
void KernelDensity<D,KT,FT,AT>::initialize_attributes(
    std::vector<DataPointType> &pts) {

  // normalize point weights
  normalize_weights(pts);

  // set masses to equal normalized weights
  for (auto &p : pts) { 
    p.attributes().set_mass(p.attributes().weight());
  }
}

template<int D, typename KT, typename FT, typename AT>
void KernelDensity<D,KT,FT,AT>::normalize_weights(std::vector<DataPointType> &pts) {

  FloatType weight_total = ConstantTraits<FloatType>::zero();
  for (const auto &p : pts) { weight_total += p.attributes().weight(); }
  for (auto &p : pts) { 
    p.attributes().set_weight(
        p.attributes().weight() / weight_total
    );
  }

}

template<int D, typename KT, typename FT, typename AT>
KernelDensity<D,KT,FT,AT>::~KernelDensity() {}



// user wrapper for single tree kde evaluation. 
// computes with the default kernel. 
template<int D, typename KT, typename FT, typename AT>
typename KernelDensity<D,KT,FT,AT>::FloatType
KernelDensity<D,KT,FT,AT>::eval(DataPointType &p, 
    FloatType rel_err, FloatType abs_err) const {

  FloatType result = eval(p.point(), kernel_, rel_err, abs_err);

  p.attributes().set_upper(result); 
  p.attributes().set_lower(result);

  return result;

}

// single point kde evaluation. based on the following algorithms:
// + ''Multiresolution Instance-Based Learning'' by Deng and Moore
// + ''Nonparametric Density Estimation: Toward Computational Tractability'' 
//   by Gray and Moore
template<int D, typename KT, typename FT, typename AT>
  template<typename KernT>
typename KernelDensity<D,KT,FT,AT>::FloatType
KernelDensity<D,KT,FT,AT>::eval(
    const GeomPointType &p, 
    const KernT &kernel,
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
  FloatType normalization = kernel.normalization();

  single_tree(data_tree_.root_, p, kernel,
              upper, lower, du, dl, 
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


template<int D, typename KT, typename FT, typename AT>
  template<typename KernT>
void KernelDensity<D,KT,FT,AT>::single_tree(
    const TreeNodeType *D_node, const GeomPointType &p, const KernT &kernel,
    FloatType &upper, FloatType &lower, 
    FloatType du, FloatType dl, 
    FloatType rel_err, FloatType abs_err) const {

  // update the kernel contributions due to points in `D_node` 
  // towards the upper/lower bounds on the kde value at point `p`. 
  FloatType du_new, dl_new; 
  estimate_contributions(D_node, p, kernel, du_new, dl_new);

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

    single_tree_base(D_node, p, kernel, du, dl, upper, lower);

  // branch: case 2: non-leaf. recursively tighten the bounds. 
  } else {

    // tighten the bounds for faster convergence
    tighten_bounds(D_node, du_new, dl_new, du, dl, upper, lower);

    // decide which halfspace is closer to the query
    TreeNodeType *closer = D_node->left_, *further = D_node->right_;
    apply_closer_heuristic(&closer, &further, p);
    
    // recursively tighten the bounds, closer halfspace first 
    single_tree(closer, p, kernel, upper, lower, du_new, dl_new, rel_err, abs_err);
    single_tree(further, p, kernel, upper, lower, du_new, dl_new, rel_err, abs_err);

  }
}

// input invariants:
// + lower <= upper, dl <= du
//
// output invariants:
// + lower <= upper
template<int D, typename KT, typename FT, typename AT>
  template<typename KernT>
void KernelDensity<D,KT,FT,AT>::single_tree_base(
    const TreeNodeType *D_node, const GeomPointType &p, const KernT &kernel,
    FloatType du, FloatType dl, 
    FloatType &upper, FloatType &lower) const {

  FloatType delta;
  for (auto i = D_node->start_idx_; i <= D_node->end_idx_; ++i) {

    delta = kernel.unnormalized_eval(p, data_tree_.points_[i].point(), 
                                        data_tree_.points_[i].attributes().lower_abw());

    delta *= data_tree_.points_[i].attributes().weight();
    upper += delta; lower += delta;
  }
  upper -= D_node->attr_.weight() * du; 
  lower -= D_node->attr_.weight() * dl;

  // see comment in tighten_bounds. 
  if (lower > upper) { upper = lower; }

}


// user wrapper for tree multi-point kernel density evaluation.
// computes with the default kernel. 
template<int D, typename KT, typename FT, typename AT>
void KernelDensity<D,KT,FT,AT>::eval(

#ifndef __CUDACC__
    std::vector<DataPointType> &queries, 
    FloatType rel_err, FloatType abs_err, 
    int leaf_nmax
#else
    std::vector<DataPointType> &queries, 
    FloatType rel_err, FloatType abs_err, 
    int leaf_nmax, 
    size_t block_size
#endif
    
    ) const {


  // construct a query tree
  KdtreeType query_tree(std::move(queries), leaf_nmax);

#ifndef __CUDACC__
  eval(query_tree, kernel_, rel_err, abs_err);
#else
  eval(query_tree, kernel_, rel_err, abs_err, block_size);
#endif

  // move the results back
  queries = std::move(query_tree.points_);

}

// user wrapper for tree multi-point kernel density evaluation.
// computes with the default kernel. 
template<int D, typename KT, typename FT, typename AT>
inline void KernelDensity<D,KT,FT,AT>::eval(

#ifndef __CUDACC__
    KdtreeType &query_tree, 
    FloatType rel_err, FloatType abs_err
#else
    KdtreeType &query_tree,
    FloatType rel_err, FloatType abs_err, 
    size_t block_size
#endif
    
    ) const {

#ifndef __CUDACC__
  eval(query_tree, kernel_, rel_err, abs_err);
#else
  eval(query_tree, kernel_, rel_err, abs_err, block_size);
#endif

}


// tree multi-point kde evaluation. computes with arbitrary kernels.
template<int D, typename KT, typename FT, typename AT>
  template<typename KernT>
void KernelDensity<D,KT,FT,AT>::eval(

#ifndef __CUDACC__
    KdtreeType &query_tree, const KernT &kernel,
    FloatType rel_err, FloatType abs_err
#else
    KdtreeType &query_tree, const KernT &kernel,
    FloatType rel_err, FloatType abs_err, 
    size_t block_size
#endif
    
    ) const {

  // initialize upper/lower bounds of individual queries to be
  // such that all data points contribute maximally/minimally
  for (auto &q : query_tree.points_) { 
    q.attributes().set_lower(0);
    q.attributes().set_upper(data_tree_.root_->attr_.weight());
  }
  query_tree.refresh_node_attributes(query_tree.root_);

  FloatType du = 1.0, dl = 0.0;


  // dual tree algorithm
  FloatType normalization = kernel.normalization(); 

#ifndef __CUDACC__
  dual_tree(data_tree_.root_, query_tree.root_, kernel,
            du, dl, rel_err, abs_err/normalization, query_tree);
#else
  CudaDirectKde<D,KernelFloatType,KernT> 
    cu_kde(data_tree_.points(), query_tree.points());
  cu_kde.kernel() = kernel;

  std::vector<KernelFloatType> host_result_cache(query_tree.size());

  dual_tree(data_tree_.root_, query_tree.root_, kernel,
            du, dl, rel_err, abs_err/normalization, query_tree,
            cu_kde, host_result_cache, block_size);
#endif

  // remember to normalize
  for (auto &q : query_tree.points_) { 

    q.attributes().set_lower(q.attributes().lower()*normalization);
    q.attributes().set_upper(q.attributes().upper()*normalization);

    report_error(std::cerr, q.point(), 
                 q.attributes().upper(), q.attributes().lower(), 
                 rel_err, abs_err);
  }

  return;
}


// tighten the contribution from all points in D_node to the upper/lower
// bounds of Q_node as well as each individual queries in Q_node
// `du`, `dl` are the present kernel contributions of D_node to bounds in Q_node. 
//
// the lower/upper bounds of Q_node is the min/max of all lower/upper 
// bounds of the individual queries 
template<int D, typename KT, typename FT, typename AT>
  template<typename KernT>
void KernelDensity<D,KT,FT,AT>::dual_tree(

#ifndef __CUDACC__
    const TreeNodeType *D_node, TreeNodeType *Q_node, const KernT &kernel,
    FloatType du, FloatType dl, FloatType rel_err, FloatType abs_err,
    KdtreeType &query_tree
#else
    const TreeNodeType *D_node, TreeNodeType *Q_node, const KernT &kernel,
    FloatType du, FloatType dl, FloatType rel_err, FloatType abs_err,
    KdtreeType &query_tree,
    CudaDirectKde<D,KernelFloatType,KernT> &cu_kde,
    std::vector<KernelFloatType> &host_result_cache,
    size_t block_size
#endif
    
    ) const {

  // update the kernel contributions due to D_node
  FloatType du_new, dl_new;
  estimate_contributions(D_node, Q_node->bbox_, kernel, du_new, dl_new);

  // BOUND: decide whether the approximation satsifies the error guarantees
  if (can_approximate(D_node, Q_node, 
        du_new, dl_new, du, dl, rel_err, abs_err)) {

    // tighten the lower/upper bound of Q_node itself
    tighten_bounds(D_node, Q_node, du_new, dl_new, du, dl);

    // tighten the individual queries
    FloatType upper_q, lower_q;
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

#ifndef __CUDACC__
    dual_tree_base(D_node, Q_node, kernel, du, dl, query_tree);
#else
    dual_tree_base(D_node, Q_node, kernel, du, dl, query_tree, 
                   cu_kde, host_result_cache, block_size);
#endif
    
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

#ifndef __CUDACC__
      dual_tree(closer, Q_node, kernel, 
          du_new, dl_new, rel_err, abs_err, query_tree);
      dual_tree(further, Q_node, kernel, 
          du_new, dl_new, rel_err, abs_err, query_tree);
#else
      dual_tree(closer, Q_node, kernel,
          du_new, dl_new, rel_err, abs_err, query_tree,
          cu_kde, host_result_cache, block_size);
      dual_tree(further, Q_node, kernel,
          du_new, dl_new, rel_err, abs_err, query_tree,
          cu_kde, host_result_cache, block_size);
#endif

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

#ifndef __CUDACC__
        dual_tree(D_node, Q_node->left_, kernel, 
            du_new, dl_new, rel_err, abs_err, query_tree);
        dual_tree(D_node, Q_node->right_, kernel, 
            du_new, dl_new, rel_err, abs_err, query_tree);
#else 
        dual_tree(D_node, Q_node->left_, kernel,
            du_new, dl_new, rel_err, abs_err, query_tree,
            cu_kde, host_result_cache, block_size);
        dual_tree(D_node, Q_node->right_, kernel,
            du_new, dl_new, rel_err, abs_err, query_tree,
            cu_kde, host_result_cache, block_size);
#endif

      // case 3: neither Q nor D are leaves
      } else {

        // tighten Q->left
        TreeNodeType *closer = D_node->left_, *further = D_node->right_;
        apply_closer_heuristic(&closer, &further, Q_node->left_->bbox_);

#ifndef __CUDACC__
        dual_tree(closer, Q_node->left_, kernel, 
            du_new, dl_new, rel_err, abs_err, query_tree);
        dual_tree(further, Q_node->left_, kernel, 
            du_new, dl_new, rel_err, abs_err, query_tree);
#else
        dual_tree(closer, Q_node->left_, kernel,
            du_new, dl_new, rel_err, abs_err, query_tree,
            cu_kde, host_result_cache, block_size);
        dual_tree(further, Q_node->left_, kernel,
            du_new, dl_new, rel_err, abs_err, query_tree,
            cu_kde, host_result_cache, block_size);
#endif

        // tighten Q->right
        closer = D_node->left_; further = D_node->right_;
        apply_closer_heuristic(&closer, &further, Q_node->right_->bbox_);

#ifndef __CUDACC__
        dual_tree(closer, Q_node->right_, kernel, 
            du_new, dl_new, rel_err, abs_err, query_tree);
        dual_tree(further, Q_node->right_, kernel, 
            du_new, dl_new, rel_err, abs_err, query_tree);
#else
        dual_tree(closer, Q_node->right_, kernel,
            du_new, dl_new, rel_err, abs_err, query_tree,
            cu_kde, host_result_cache, block_size);
        dual_tree(further, Q_node->right_, kernel,
            du_new, dl_new, rel_err, abs_err, query_tree,
            cu_kde, host_result_cache, block_size);
#endif

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


template<int D, typename KT, typename FT, typename AT>
  template<typename KernT>
void KernelDensity<D,KT,FT,AT>::dual_tree_base(
#ifndef __CUDACC__
    const TreeNodeType *D_node, TreeNodeType *Q_node, const KernT &kernel,
    FloatType du, FloatType dl, 
    KdtreeType &query_tree
#else
    const TreeNodeType *D_node, TreeNodeType *Q_node, const KernT &kernel,
    FloatType du, FloatType dl, 
    KdtreeType &query_tree,
    CudaDirectKde<D,KernelFloatType,KernT> &cu_kde,
    std::vector<KernelFloatType> &host_result_cache,
    size_t block_size
#endif
    ) const {

#ifdef __CUDACC__
  cu_kde.unnormalized_eval(
      D_node->start_idx_, D_node->end_idx_, 
      Q_node->start_idx_, Q_node->end_idx_, 
      host_result_cache, block_size);
#endif

  FloatType min_q = std::numeric_limits<FloatType>::max();
  FloatType max_q = std::numeric_limits<FloatType>::min();

  FloatType lower_q, upper_q;
  for (auto i = Q_node->start_idx_; i <= Q_node->end_idx_; ++i) {

    // update the contribution of each point due to D_node
    upper_q = query_tree.points_[i].attributes().upper();
    lower_q = query_tree.points_[i].attributes().lower();

#ifndef __CUDACC__

    single_tree_base(
        D_node, query_tree.points_[i].point(), kernel,
        1.0, 0.0, upper_q, lower_q);

#else

    upper_q += host_result_cache[i-(Q_node->start_idx_)]; 
    upper_q -= D_node->attr_.weight();

    lower_q += host_result_cache[i-(Q_node->start_idx_)]; 

    if (lower_q > upper_q) { upper_q = lower_q; }

#endif

    query_tree.points_[i].attributes().set_lower(lower_q);
    query_tree.points_[i].attributes().set_upper(upper_q);

    min_q = std::min(lower_q, min_q);
    max_q = std::max(upper_q, max_q);

  }

  Q_node->attr_.set_lower(min_q);
  Q_node->attr_.set_upper(max_q);

}


template<int D, typename KT, typename FT, typename AT>
void KernelDensity<D,KT,FT,AT>::tighten_bounds(
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
template<int D, typename KT, typename FT, typename AT>
void KernelDensity<D,KT,FT,AT>::tighten_bounds(
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



template<int D, typename KT, typename FT, typename AT>
inline bool KernelDensity<D,KT,FT,AT>::can_approximate(
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
template<int D, typename KT, typename FT, typename AT>
bool KernelDensity<D,KT,FT,AT>::can_approximate(
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




template<int D, typename KT, typename FT, typename AT>
  template<typename ObjT>
inline void KernelDensity<D,KT,FT,AT>::apply_closer_heuristic(
    TreeNodeType **closer, TreeNodeType **further, const ObjT &obj) const {

  if ((*closer)->bbox_.min_dist(obj) > (*further)->bbox_.min_dist(obj)) {
    std::swap(*closer, *further);
  }

}


template<int D, typename KT, typename FT, typename AT>
  template<typename ObjT, typename KernT> 
void KernelDensity<D,KT,FT,AT>::estimate_contributions(
    const TreeNodeType *D_node, const ObjT &obj, const KernT &kernel,
    FloatType &du, FloatType &dl) const {

  GeomPointType proxy;
  const static GeomPointType origin;

  // use the minimum(maximum) distance to the argument in each 
  // dimension to bound the min/max kernel contributions

  for (int i = 0; i < D; ++i) { proxy[i] = D_node->bbox_.min_dist(i, obj); }
  du = kernel.unnormalized_eval(proxy, origin, D_node->attr_.upper_abw());

  for (int i = 0; i < D; ++i) { proxy[i] = D_node->bbox_.max_dist(i, obj); }
  dl = kernel.unnormalized_eval(proxy, origin, D_node->attr_.lower_abw());

}


template<int D, typename KT, typename FT, typename AT>
void KernelDensity<D,KT,FT,AT>::report_error(
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

// user wrapper for direct kernel density evaluation.
// computes with the default kernel. 
template<int D, typename KT, typename FT, typename AT>
typename KernelDensity<D,KT,FT,AT>::FloatType
KernelDensity<D,KT,FT,AT>::direct_eval(DataPointType &p) const {
  FloatType result = direct_eval(p.point(), kernel_);
  p.attributes().set_upper(result);
  p.attributes().set_lower(result);
  return result;
}

// direct kernel density evaluation. computes using arbitrary kernels. 
template<int D, typename KT, typename FT, typename AT>
  template<typename KernT>
typename KernelDensity<D,KT,FT,AT>::FloatType
KernelDensity<D,KT,FT,AT>::direct_eval(
    const GeomPointType &p, const KernT &kernel) const {

  FloatType total = ConstantTraits<FloatType>::zero();
  for (const auto &datum : data_tree_.points()) {
    total += 
      datum.attributes().weight() * 
      kernel.unnormalized_eval(p, datum.point(), 
                                  datum.attributes().lower_abw() );
  }
  total *= kernel.normalization();
  return total;

}


// user wrapper for direct kernel density evaluation.
#ifndef __CUDACC__
template<int D, typename KT, typename FT, typename AT>
inline void KernelDensity<D,KT,FT,AT>::direct_eval(
    std::vector<DataPointType> &queries) const {
  direct_eval(queries, kernel_);
  return; 
}
#else 
template<int D, typename KT, typename FT, typename AT>
inline void KernelDensity<D,KT,FT,AT>::direct_eval(
    std::vector<DataPointType> &queries, size_t block_size
    ) const {
  direct_eval(queries, kernel_, block_size);
  return; 
}
#endif



// direct kernel density evaluation. computes using arbitrary kernels. 
template<int D, typename KT, typename FT, typename AT>
  template<typename KernT>
void KernelDensity<D,KT,FT,AT>::direct_eval(

#ifndef __CUDACC__
    std::vector<DataPointType> &queries, const KernT &kernel
#else
    std::vector<DataPointType> &queries, const KernT &kernel, size_t block_size
#endif

    ) const {

#ifndef __CUDACC__
  for (auto &q : queries) {
    FloatType result = direct_eval(q.point(), kernel);
    q.attributes().set_lower(result);
    q.attributes().set_upper(result);
  }
#else
  std::vector<KernelFloatType> host_results(queries.size());

  CudaDirectKde<D,KernelFloatType,KernT> cuda_kde(data_tree_.points(), queries);
  cuda_kde.kernel() = kernel;

  cuda_kde.eval(0, data_tree_.size()-1, 0, queries.size()-1, 
                host_results, block_size);

  for (size_t i = 0; i < queries.size(); ++i) {
    queries[i].attributes().set_lower(host_results[i]);
    queries[i].attributes().set_upper(host_results[i]);
  }
#endif

  return; 
}


}

#endif
