#ifndef BBRCITKDE_KDTREE_H__
#define BBRCITKDE_KDTREE_H__

#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include <utility>
#include <stack>
#include <iostream>

#include <DecoratedPoint.h>
#include <Rectangle.h>
#include <FloatUtils.h>
#include <Attributes/PointWeights.h>

// API
// ---

namespace bbrcit {

template<int D, typename AttrT, typename FloatT> class Kdtree;

template<int D, typename KernelT, typename FloatT, typename AttrT> class KernelDensity;

template<int D, typename AttrT, typename FloatT>
void swap(Kdtree<D,AttrT,FloatT>&, Kdtree<D,AttrT,FloatT>&);

// Kdtree<> implements a D-dimensional kdtree. It currently supports:
// + Range search. 
template<int D, 
         typename AttrT=PointWeights<int>, 
         typename FloatT = double>
class Kdtree {

  public: 

    using DataPointType = DecoratedPoint<D,AttrT,FloatT>;
    using RectangleType = Rectangle<D,FloatT>;
    using AttributesType = AttrT;
    using FloatType = FloatT;
    static constexpr int dim() { return D; }

  protected:
    using IndexType = typename std::vector<DataPointType>::size_type;
    friend void swap<>(Kdtree<D,AttrT,FloatT>&, Kdtree<D,AttrT,FloatT>&);

    template <int DIM, typename KT, typename FT, typename AT> 
      friend class KernelDensity;

  public: 

    // default constructor yields a null tree. 
    Kdtree();

    // construct Kdtree out of data given as a list of DataPointType's. 
    Kdtree(const std::vector<DataPointType> &data, int leaf_nmax=2);
    Kdtree(std::vector<DataPointType> &&data, int leaf_nmax=2);

    // copy-control. operator='s use copy and swap.
    Kdtree(const Kdtree<D,AttrT,FloatT>&);
    Kdtree(Kdtree<D,AttrT,FloatT>&&) noexcept;
    Kdtree<D,AttrT,FloatT>& operator=(Kdtree<D,AttrT,FloatT>);
    virtual ~Kdtree();

    // returns true if this is a null tree. 
    bool empty() const;

    // returns the number of data points in this Kdtree.
    IndexType size() const;

    // returns the maximum number of points per leaf. 
    int leaf_nmax() const;

    // returns a const reference to the points.
    const std::vector<DataPointType>& points() const;

    // prints each DataPointType stored in the Kdtree on a separate line os.
    void print_points(std::ostream &os) const;

    // (1) print each partition at depth d and leaf partitions of depth <= d on a separate line of os.
    // (2) returns each partition at depth d and leaf partitions of depth <= d in a vector.
    void print_partitions(int d, std::ostream&) const;
    void report_partitions(int d, std::vector<RectangleType>&) const;

    // (1) print each DataPointType contained in the query window on a separate line of os.
    // (2) returns each DataPointType contained in the query window in a vector. 
    void print_range_search(const RectangleType &query, std::ostream &os) const;
    void range_search(const RectangleType &query, std::vector<DataPointType>&) const;

    // primarily for debugging: 
    // + report_leaves: save ranges of point indices for every leaf. 
    // + root_attributes: returns the attributes object of the root node. 
    void report_leaves(std::vector<std::pair<IndexType,IndexType>>&) const;

  protected:

    // Kdtree<>::Node represents a node in the Kdtree. 
    // (I/L) are members that are meaningful for internal/leaf nodes. 
    // + An object represents an internal node iff left_=right_=nullptr.
    struct Node {

      // (I) coordinate at which to partition half spaces.
      FloatType split_ = FloatType();

      // (I/L) links to the daughters
      Node *left_ = nullptr, *right_ = nullptr;

      // (I/L) the index range [start_idx_, end_idx_] are the data points 
      // in points_ that are organized under this node
      IndexType start_idx_ = 0;
      IndexType end_idx_ = 0;

      // (I/L) attributes associated with this node
      AttributesType attr_ = AttributesType();

      // (I/L) the smallest rectangle containing all data points 
      // associated with this node. 
      RectangleType bbox_;
      
      // returns true if this object is a leaf
      bool is_leaf() const { return left_ == nullptr && right_ == nullptr; }
      
      // returns the number of points under this node 
      IndexType size() const { return end_idx_ - start_idx_ + 1;}

    };

    // list of DataPointType's that stores the data from which to construct the Kdtree. 
    // Note: the order of DataPointType's should NOT change once the object is constructed.
    std::vector<DataPointType> points_;

    // root node of the Kdtree. 
    Node *root_;

    // maximum number of points per leaf. 
    int leaf_nmax_;

  private:

    // helper functions
    void initialize();
    void merge_duplicates();
    RectangleType compute_bounding_box() const;
    Node* construct_tree(int, int, int, const RectangleType&);

    Node* tree_deep_copy(const Node*) const;

    void delete_tree(Node*);

    void retrieve_point_indices(const Node*, std::vector<IndexType>&) const;
    void retrieve_range_indices(const Node*, const RectangleType&, std::vector<IndexType>&) const;
    void retrieve_partitions(const Node *, int, std::vector<RectangleType>&) const;

    void refresh_node_attributes(Node*);
};

// Implementations
// ---------------

// refresh node attributes stored in the subtree pointed to by `p`, 
// according to the current attributes of its data points (`points`). 
template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::refresh_node_attributes(Node *p) {

  if (p == nullptr) { return; }

  if (p->is_leaf()) {
    size_t i = p->start_idx_, j = p->end_idx_;
    p->attr_ = points_[i].attributes();
    for (int k = i+1; k <= j; ++k) {
      p->attr_.merge(points_[k].attributes());
    }
  } else {
    refresh_node_attributes(p->left_);
    refresh_node_attributes(p->right_);
    p->attr_ = merge(p->left_->attr_, p->right_->attr_);
  }

  return;
}

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::print_partitions(int depth, std::ostream &os) const {
  std::vector<RectangleType> result;
  retrieve_partitions(root_, depth, result);
  for (const auto &r : result) { os << r << std::endl; }
}

template<int D, typename AttrT, typename FloatT>
inline void Kdtree<D,AttrT,FloatT>::report_partitions(int depth, std::vector<RectangleType> &result) const {
  retrieve_partitions(root_, depth, result);
}

// use DFS to retrieve partitions. 
template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::retrieve_partitions(
    const Node *r, int depth, 
    std::vector<RectangleType> &result) const {

  // handle null tree separately
  if (r == nullptr) { return; }

  // both conditions required. see API. 
  if (depth == 0 || r->is_leaf()) { result.push_back(r->bbox_); return; }

  retrieve_partitions(r->left_, depth-1, result);
  retrieve_partitions(r->right_, depth-1, result);

}

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::print_range_search(const RectangleType &query_range, std::ostream &os) const {

  std::vector<IndexType> result_indices;
  retrieve_range_indices(root_, query_range, result_indices);
  for (auto i : result_indices) { os << points_[i] << std::endl; }

}

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::range_search(const RectangleType &query_range, std::vector<DataPointType> &result) const {
  std::vector<IndexType> result_indices;
  retrieve_range_indices(root_, query_range, result_indices);
  for (auto i : result_indices) { result.emplace_back(points_[i]); }
}

// standard range search algorithm for kdtrees. 
template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::retrieve_range_indices(
    const Node *v, const RectangleType &query_range, 
    std::vector<IndexType> &result) const {

  if (v == nullptr) { return; }

  if (v->is_leaf()) { 
    for (IndexType i = v->start_idx_; i <= v->end_idx_; ++i) {
      if (query_range.contains(points_[i])) {
        result.push_back(i);
      }
    }
  } else {

    // left halfspace
    if (query_range.contains(v->left_->bbox_)) {
      retrieve_point_indices(v->left_, result);
    } else {
      retrieve_range_indices(v->left_, query_range, result);
    }

    // right halfspace
    if (query_range.contains(v->right_->bbox_)) {
      retrieve_point_indices(v->right_, result);
    } else {
      retrieve_range_indices(v->right_, query_range, result);
    }

  }
}

template<int D, typename AttrT, typename FloatT>
void swap(Kdtree<D,AttrT,FloatT> &lhs, Kdtree<D,AttrT,FloatT> &rhs) {

  using std::swap;
  
  swap(lhs.points_, rhs.points_);
  swap(lhs.leaf_nmax_, rhs.leaf_nmax_);

  // simply switch the root_ pointers.
  typename Kdtree<D,AttrT,FloatT>::Node *ptemp = lhs.root_;
  lhs.root_ = rhs.root_;
  rhs.root_ = ptemp;
}

template<int D, typename AttrT, typename FloatT>
inline Kdtree<D,AttrT,FloatT>& Kdtree<D,AttrT,FloatT>::operator=(Kdtree<D,AttrT,FloatT> rhs) {
  swap(*this, rhs); return *this;
}

template<int D, typename AttrT, typename FloatT>
inline bool Kdtree<D,AttrT,FloatT>::empty() const { return root_ == nullptr; }

template<int D, typename AttrT, typename FloatT>
inline typename Kdtree<D,AttrT,FloatT>::IndexType Kdtree<D,AttrT,FloatT>::size() const { 
  return root_ == nullptr ? 0 : root_->size(); 
}

template<int D, typename AttrT, typename FloatT>
inline int Kdtree<D,AttrT,FloatT>::leaf_nmax() const 
{ return leaf_nmax_; }

// DFS to the leaves and print the index range of points to os
template<int D, typename AttrT, typename FloatT>
inline void Kdtree<D,AttrT,FloatT>::report_leaves(
    std::vector<std::pair<IndexType,IndexType>> &result) const { 
  std::stack<Node*> s; s.push(root_);
  while (!s.empty()) {
    Node *r = s.top(); s.pop();
    if (r->is_leaf()) { 
      result.push_back(std::make_pair(r->start_idx_, r->end_idx_));
    } else {
      if (r->right_) { s.push(r->right_); }
      if (r->left_) { s.push(r->left_); }
    }
  }
}

template<int D, typename AttrT, typename FloatT>
inline const std::vector<typename Kdtree<D,AttrT,FloatT>::DataPointType>&
Kdtree<D,AttrT,FloatT>::points() const {
  return points_;
}

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::print_points(std::ostream &os) const {
  for (const auto &p : points_) { os << p << std::endl; }
}

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::retrieve_point_indices(const Node *r, std::vector<IndexType> &result) const {

  // check for null tree
  if (r == nullptr) { return; }

  // simply append all points under this node
  for (IndexType i = r->start_idx_; i <= r->end_idx_; ++i) {
    result.emplace_back(i);
  }
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree() : points_(0), root_(nullptr) {}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::~Kdtree() { delete_tree(root_); }

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::initialize() {

  // merge duplicate keys
  merge_duplicates();

  // build the tree
  if (!points_.empty()) { 
    root_ = construct_tree(0, points_.size()-1, 0, compute_bounding_box()); 
  }
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree(
    const std::vector<DataPointType> &points, int leaf_nmax) 
  : points_(points), root_(nullptr), leaf_nmax_(leaf_nmax) {
  initialize();
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree(
    std::vector<DataPointType> &&points, int leaf_nmax) 
  : points_(std::move(points)), root_(nullptr), leaf_nmax_(leaf_nmax) {
  initialize();
}

// merge duplicate keys in the data by merging the attributes. 
template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::merge_duplicates() {

  if (points_.empty()) { return; }

  // preprocess by sorting the data points lexicographically
  // note: exact comparison of floating point is ok for this purpose
  std::sort(points_.begin(), points_.end(), ExactLexicoLess<DataPointType>);

  // remove duplicates by merging attributes. 
  // algorithm is similar to the partition step in quicksort.
  size_t i = 0;
  for (size_t j = 1; j < points_.size(); ++j) {
    if (!ExactEqual(points_[i], points_[j])) {
      swap(points_[++i], points_[j]);
    } else {
      points_[i].set_attributes(
          merge(points_[i].attributes(), 
                points_[j].attributes())
      );
    }
  }

  // with duplicates removed, ok to shrink data size
  points_.resize(i+1);

}


// if points_ is non-empty, return the minimum area axis-aligned rectangle that 
// containing all DataPointType's in points_. otherwise return a degenerate rectangle 
// (a rectangle of 0 length edges at the origin)
template<int D, typename AttrT, typename FloatT>
typename Kdtree<D,AttrT,FloatT>::RectangleType Kdtree<D,AttrT,FloatT>::compute_bounding_box() const {

  // handle the empty data set. 
  if (points_.empty()) { return RectangleType(); }

  // simply scan through points_ to find the mininimum and maximum 
  // coordinate value in each dimension. 
  std::vector<FloatType> min_coord_val(D, std::numeric_limits<FloatType>::max());
  std::vector<FloatType> max_coord_val(D, std::numeric_limits<FloatType>::min());
  for (const auto &p : points_) {
    for (int i = 0; i < D; ++i) { 
      min_coord_val[i] = std::min(p[i], min_coord_val[i]);
      max_coord_val[i] = std::max(p[i], max_coord_val[i]);
    }
  }

  // initialize and resize a rectangle to the min/max coordinates. 
  RectangleType r;
  for (int i = 0; i < D; ++i) { 
    r.resize(i, {min_coord_val[i], max_coord_val[i]});
  }

  return r;
}

// returns a pointer to the root of a Kdtree constructed over elements in points_
// within the *closed* indices interval [i,j]. Uses a pre-order tree walk. 
// + d indexes the dimension to perform the first split.
// + bbox is a minimum area rectangle containing all the points in [i,j].
// + this procedure rearranges points_ such that the indices stored at the leaves
//   refer to the appropriate data point. 
// Note: do NOT change the ordering of points_ outside of this function. 
template<int D, typename AttrT, typename FloatT>
typename Kdtree<D,AttrT,FloatT>::Node* 
Kdtree<D,AttrT,FloatT>::construct_tree(
    int i, int j, int d, const RectangleType &bbox) {

  Node *p = new Node();
  p->start_idx_ = i;
  p->end_idx_ = j;

  // explicitly store the bounding box at the node. 
  p->bbox_ = bbox;

  // create a leaf node when [i,j] contains no more than leaf_nmax_ indices. 
  if (j-i+1 <= leaf_nmax_) { 

    p->attr_ = points_[i].attributes();
    for (int k = i+1; k <= j; ++k) {
      p->attr_.merge(points_[k].attributes());
    }

  // partition points_ by the lower median m. [i,m] go in the left subtree 
  // while [m+1,j] go in the right subtree. 
  } else {

    // median finding in expected linear time. Note: C++ requires the end 
    // iterator to be one past the last index, hence the +1 in the 3rd 
    // argument. 
    int m = i + (j-i) / 2;
    std::nth_element(points_.begin()+i, points_.begin()+m, points_.begin()+j+1, 
                    [d] (const DataPointType &p1, const DataPointType &p2) { return p1[d] < p2[d]; });
    p->split_ = points_[m][d];

    // pre-order tree walk, but first partition the bounding box. 
    RectangleType left_bbox = bbox.lower_halfspace(d, p->split_);
    p->left_ = construct_tree(i, m, (d+1)%D, left_bbox);

    RectangleType right_bbox = bbox.upper_halfspace(d, p->split_);
    p->right_ = construct_tree(m+1, j, (d+1)%D, right_bbox);

    p->attr_ = merge(p->left_->attr_, p->right_->attr_);
  }
  return p;
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree(Kdtree<D,AttrT,FloatT> &&obj) noexcept : 
  points_(std::move(obj.points_)), 
  leaf_nmax_(std::move(obj.leaf_nmax_)) {
  root_ = obj.root_;
  obj.root_ = nullptr;
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree(const Kdtree<D,AttrT,FloatT> &obj) 
  : points_(obj.points_), leaf_nmax_(obj.leaf_nmax_) {
  root_ = tree_deep_copy(obj.root_);
}

// return a pointer to a deep copy of the tree pointed to by the argument. 
template<int D, typename AttrT, typename FloatT>
typename Kdtree<D,AttrT,FloatT>::Node* 
Kdtree<D,AttrT,FloatT>::tree_deep_copy(
    const Kdtree<D,AttrT,FloatT>::Node *target_r) const {
  if (target_r == nullptr) { return nullptr; }
  Node *r = new Node(*target_r);
  r->left_ = tree_deep_copy(target_r->left_);
  r->right_ = tree_deep_copy(target_r->right_);
  return r;
}

// delete Kdtree<>::Nodes by a post order tree walk. 
template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::delete_tree(Node *p) {
  if (p == nullptr) { return; }
  if (!p->is_leaf()) { 
    delete_tree(p->left_);
    delete_tree(p->right_);
  }
  delete p; 
}

}

#endif
