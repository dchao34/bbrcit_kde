#ifndef BBRCITKDE_KERNELDENSITY_H__
#define BBRCITKDE_KERNELDENSITY_H__

#include <iostream>

#include <Kdtree.h>
#include <Attributes/AdaKdeAttributes.h>
#include <Kernels/EpanechnikovKernel.h>

#ifdef __CUDACC__
#include <CudaDirectKde.h>
#endif

namespace bbrcit {

template<int D, typename KT, typename FT, typename AT> class KernelDensity;

// custom swap for KernelDensity<>
template<int D, typename KT, typename FT, typename AT>
void swap(KernelDensity<D,KT,FT,AT>&, KernelDensity<D,KT,FT,AT>&);


// least squares cross validation using numerical integration for 2d data. 
#ifndef __CUDACC__
template<typename KT, typename FT, typename AT>
FT lsq_numint_cross_validate(
    const KernelDensity<2,KT,FT,AT> &kde, 
    FT start_x, FT end_x, int step_x, 
    FT start_y, FT end_y, int step_y,
    FT rel_err=1e-6, FT abs_err=1e-8, int qtree_leaf_nmax=2);
#else
template<typename KT, typename FT, typename AT>
FT lsq_numint_cross_validate(
    const KernelDensity<2,KT,FT,AT> &kde, 
    FT start_x, FT end_x, int step_x, 
    FT start_y, FT end_y, int step_y,
    FT rel_err=1e-6, FT abs_err=1e-8, int qtree_leaf_nmax=1024,
    size_t block_size=128);
#endif


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

    // returns a const reference to the data tree
    const KdtreeType& data_tree() const;


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


    // convert between adaptive and non-adaptive kernel density estimates.
    //
    // `adapt_density`: makes the kernel density adaptive. degree of the adaptation
    //                  depends on 0 <= `alpha` <= 1, with 0 being equivalent to no 
    //                  adaptation at all. the local bandwidth corrections are computed
    //                  using the current bandwidth setting. 
    //     
    // `unadapt_density`: makes the density non-adaptive.
    //
    // note: adapting densitities require all pair self evaluations of the reference
    // data points and is thus computationally intensive. for large problem sizes, a 
    // GPU is highly recommended. 
#ifndef __CUDACC__ 
    void adapt_density(FloatType alpha, 
        FloatType rel_err=1e-6, FloatType abs_err=1e-6);
#else
    void adapt_density(FloatType alpha, 
        FloatType rel_err=1e-6, FloatType abs_err=1e-6, size_t block_size=128);
#endif
    void unadapt_density();


    // compute the cross validation score for the current 
    // kernel configuration. likelihood and least squares are available. 
    // see non-member `lsq_convolution_cross_validate` for another flavor of
    // cross validation. 
    //
    // note: for adaptive kernels, convolution kernel based least squares cross
    // validation is not well defined. for this reason, least square cross validation 
    // for adaptive kernels should use `lsq_numint_cross_validate` instead of 
    // `lsq_convolution_cross_validate` 
#ifndef __CUDACC__ 
    FloatType likelihood_cross_validate(
        FloatType rel_err=1e-6, FloatType abs_err=1e-6) const;
    FloatType lsq_convolution_cross_validate(
        FloatType rel_err=1e-6, FloatType abs_err=1e-6) const;
#else
    FloatType likelihood_cross_validate(
        FloatType rel_err=1e-6, FloatType abs_err=1e-6, size_t block_size=128) const;
    FloatType lsq_convolution_cross_validate(
        FloatType rel_err=1e-6, FloatType abs_err=1e-6, size_t block_size=128) const;
#endif


    // simulate a point from this kernel density. 
    // `e` is a random number engine from `std::random``. 
    // (1) the result is stored in `p` with `p[i]` corresponding to component `i`. 
    // (2) the result is returned as a DataPointType with attributes default initialized. 
    template <typename RNG> 
      void simulate(RNG&, std::vector<FloatType> &p) const;
    template <typename RNG> 
      DataPointType simulate(RNG&) const;

  private:

    // internal state
    // --------------

    // the kernel function of the estimator 
    KernelType kernel_;

    // the reference points associated with this kernel density. the points 
    // are required to have the same attributes as those in AdaKdeAttributes<>. 
    // Below are some additional invariants maintained as well as some notable
    // comments:
    //
    // + weight: these are the relative importance of each data point. 
    //           throughout the lifetime of this KernelDensity<> object, 
    //           we maintain the sum over all point weights to be 1.0. 
    //
    // + mass: while similar to weight, this is the actual contribution that
    //         the point has towards kde queries. though this is usually the 
    //         same as weight, it can be different; for example, in the adaptive
    //         kernels, this is the weight multiplied by local bandwidth corrections. 
    //
    KdtreeType data_tree_;

    // cumulative weights of each data point. used for simulation. 
    std::vector<FloatType> cum_weights_;

    // helper functions for initialization
    // ------------------------------------------
    void initialize_attributes(std::vector<DataPointType>&);
    void normalize_weights(std::vector<DataPointType>&);
    void initialize_cum_weights();


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

}

#include "KernelDensityImpl.h"

#endif
