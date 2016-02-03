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

// KernelDensity<> implements a kernel density estimator in D dimensions 
// using a D dimensional Kdtree. 
template<int D, 
         typename KernelT=EpanechnikovKernel<D>,
         typename DAttrT=KdeAttributes<double>,
         typename FloatT=double>
class KernelDensity {

  private:

    using KernelDensityType = KernelDensity<D,KernelT,DAttrT,FloatT>;

    using DataTreeType = Kdtree<D,DAttrT,FloatT>; 
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

    // return the kde evaluated at point `p`. its relative error is at most `rel_err`.
    FloatT eval(const GeomPointType &p, FloatType rel_err) const;

    // return the kde evaluate at points in `queries`. the relative error is at most `rel_err`.
    void eval(std::vector<DataPointType> &queries, FloatType rel_err) const;


  private:

    // bandwidth of the estimator
    FloatType bandwidth_;

    // the kernel function of the estimator 
    KernelType kernel_;

    // the kernel function of the estimator 
    DataTreeType data_tree_;


    // single tree traversal for the single point eval()
    void single_tree(DataNodeType*, const GeomPointType&, 
                     FloatType&, FloatType&, FloatType, FloatType, 
                     FloatType, 
                     const GeomRectangleType&, size_t) const;

    // dual tree traversal for the multi point eval()
    void dual_tree(DataNodeType*, DataNodeType*, 
                    FloatType, FloatType, FloatType, 
                    const GeomRectangleType&, size_t, 
                    const GeomRectangleType&, size_t, 
                    DataTreeType&) const;

  // Mainly for debugging: 
  // ---------------------
  // + naive kde evaluation; slow... O(n^2). 
  // + eval_iterative: evaluate the kde at point `p`. 
  //                   relative error is no more than `rel_err`,
  //                   relative cutoff tolerance is `rel_tol`.
  public:

    FloatT naive_eval(const GeomPointType&) const;
    FloatT eval_iterative(const GeomPointType &p, FloatType rel_err, FloatType rel_tol) const;

  private:

    // priority queue node used in eval_iterative(). 
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


// single point kde evaluation
template<int D, typename KT, typename DT, typename FT>
typename KernelDensity<D,KT,DT,FT>::FloatType
KernelDensity<D,KT,DT,FT>::eval(const GeomPointType &p, FloatType rel_err) const {

  FloatType upper = data_tree_.root_->attr_.weight();
  FloatType lower = ConstantTraits<FloatType>::zero();
  FloatType du = 1.0, dl = 0.0;

  single_tree(data_tree_.root_, p, 
              upper, lower, du, dl, rel_err, 
              data_tree_.bbox_, 0);

  FloatType result = lower + (upper - lower) / 2;
  result *= KernelType::normalization;
  result /= (std::pow(bandwidth_, D) * data_tree_.size());
  return result;

}

template<int D, typename KT, typename DT, typename FT>
void KernelDensity<D,KT,DT,FT>::single_tree(
    DataNodeType *r, const GeomPointType &p, 
    FloatType &upper, FloatType &lower, FloatType du, FloatType dl, 
    FloatType rel_err, 
    const GeomRectangleType &bbox, size_t depth) const {

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
  //if (approxmately_equal(du_new, dl_new, rel_err, rel_err)) { return; }
  if (std::abs(du_new - dl_new) < 2 * (lower + dl_new) * rel_err / data_tree_.size()) {
    return; 
  }

  // exclusion pruning
  if (almost_equal(du_new, ConstantTraits<FloatType>::zero())) { return; }

  if (r->is_leaf()) {
    for (auto i = r->start_idx_; i <= r->end_idx_; ++i) {
      FloatType delta = kernel_.unnormalized_eval(
          (p - data_tree_.points_[i].point()) / bandwidth_
      );
      upper += delta; lower += delta;
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
    single_tree(closer, p, 
                upper, lower, du_new, dl_new, rel_err, 
                *closer_r, (depth+1)%D);
    single_tree(further, p, 
                upper, lower, du_new, dl_new, rel_err, 
                *further_r, (depth+1)%D);

  }
}


// multi-point kde evaluation
template<int D, typename KT, typename DT, typename FT>
void KernelDensity<D,KT,DT,FT>::eval(
    std::vector<DataPointType> &queries, FloatType rel_err) const {

  for (auto &q : queries) { 
    q.attributes().set_lower(0);
    q.attributes().set_upper(data_tree_.size());
  }
  FloatType du = 1.0, dl = 0.0;

  DataTreeType query_tree(std::move(queries));

  dual_tree(query_tree.root_, data_tree_.root_, 
            du, dl, rel_err,
            query_tree.bbox_, 0, 
            data_tree_.bbox_, 0, 
            query_tree);

  FloatType normalization = KernelType::normalization; 
  normalization /= (std::pow(bandwidth_, D) * data_tree_.size());
  for (auto &q : query_tree.points_) { 
    q.attributes().set_lower(q.attributes().lower()*normalization);
    q.attributes().set_upper(q.attributes().upper()*normalization);
  }

  queries = std::move(query_tree.points_);

  return;

}


template<int D, typename KT, typename DT, typename FT>
void KernelDensity<D,KT,DT,FT>::dual_tree(
    DataNodeType *Q_node, DataNodeType *D_node, 
    FloatType du, FloatType dl, FloatType rel_tol, 
    const GeomRectangleType &Q_box, size_t Q_depth, 
    const GeomRectangleType &D_box, size_t D_depth, 
    DataTreeType &query_tree) const {

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
      auto &query_attr = query_tree.points_[i].attributes();
      query_attr.set_lower(query_attr.lower() + D_weight * dl_new);
      query_attr.set_upper(query_attr.upper() + D_weight * (du_new - 1));
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
      auto &query_attr = query_tree.points_[i].attributes();
      for (auto j = D_node->start_idx_; j <= D_node->end_idx_; ++j) {
        FloatType delta = kernel_.unnormalized_eval(
          (query_tree.points_[i].point() - data_tree_.points_[j].point()) / bandwidth_
        );
        query_attr.set_lower(query_attr.lower() + delta);
        query_attr.set_upper(query_attr.upper() + delta);
      }
      query_attr.set_upper(query_attr.upper() - D_weight);

      min_q = std::min(query_attr.lower(), min_q);
      max_q = std::max(query_attr.upper(), max_q);

    }

    Q_node->attr_.set_lower(min_q);
    Q_node->attr_.set_upper(max_q);

  } else {

    // case 2: Q is leaf
    if (Q_node->is_leaf()) {

      GeomRectangleType D_lbox = D_box.lower_halfspace(D_depth, D_node->split_);
      GeomRectangleType D_ubox = D_box.upper_halfspace(D_depth, D_node->split_);

      // closer heuristic
      DataNodeType *closer = D_node->left_; const GeomRectangleType *closer_r = &D_lbox;
      DataNodeType *further = D_node->right_; const GeomRectangleType *further_r = &D_ubox;
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
        DataNodeType *closer = D_node->left_; const GeomRectangleType *closer_r = &D_lbox;
        DataNodeType *further = D_node->right_; const GeomRectangleType *further_r = &D_ubox;
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





// ---------------------------------------------------------------
// ---------------------------------------------------------------

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
typename KernelDensity<D,KT,DT,FT>::FloatType
KernelDensity<D,KT,DT,FT>::eval_iterative(const GeomPointType &p, 
                                   FloatType rel_err, FloatType rel_tol) const {

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


}

#endif
