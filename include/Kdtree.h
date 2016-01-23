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
#include <PointWeights.h>
#include <EmptyNodeAttributes.h>

// API
// ---

namespace bbrcit {

template<int D, typename PAttrT, typename NAttrT, typename FloatT> class Kdtree;

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void swap(Kdtree<D,PAttrT,NAttrT,FloatT>&, Kdtree<D,PAttrT,NAttrT,FloatT>&);

// Kdtree<> implements a D-dimensional kdtree. It currently supports:
// + Range search. 
template<int D, 
         typename PAttrT=PointWeights<int>, 
         typename NAttrT=EmptyNodeAttributes, 
         typename FloatT = double>
class Kdtree {

  public: 

    using DataPointType = DecoratedPoint<D,PAttrT,FloatT>;
    using FloatType = FloatT;
    using RectangleType = Rectangle<D,FloatType>;
    using NodeAttributesType = NAttrT;

    static constexpr int dim() { return D; }

  protected:
    using IndexType = typename std::vector<DataPointType>::size_type;
    friend void swap<>(Kdtree<D,PAttrT,NAttrT,FloatT>&, Kdtree<D,PAttrT,NAttrT,FloatT>&);

  public: 

    // default constructor yields a null tree. 
    Kdtree();

    // construct Kdtree out of data given as a list of DataPointType's. 
    Kdtree(const std::vector<DataPointType> &data, int leaf_nmax=2);
    Kdtree(std::vector<DataPointType> &&data, int leaf_nmax=2);

    // copy-control. operator='s use copy and swap.
    Kdtree(const Kdtree<D,PAttrT,NAttrT,FloatT>&);
    Kdtree(Kdtree<D,PAttrT,NAttrT,FloatT>&&) noexcept;
    Kdtree<D,PAttrT,NAttrT,FloatT>& operator=(Kdtree<D,PAttrT,NAttrT,FloatT>);
    virtual ~Kdtree();

    // returns true if this is a null tree. 
    bool empty() const;

    // returns the number of data points in this Kdtree.
    IndexType size() const;

    // returns the maximum number of points per leaf. 
    int leaf_nmax() const;

    // returns the bounding box of the points.
    const RectangleType& bounding_box() const;

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
    const NAttrT& root_attributes() const { return root_->attr_; }

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
      NodeAttributesType attr_ = NodeAttributesType();
      
      // returns true if this object is a leaf
      bool is_leaf() const { return left_ == nullptr && right_ == nullptr; }

    };

    // list of DataPointType's that stores the data from which to construct the Kdtree. 
    // Note: the order of DataPointType's should NOT change once the object is constructed.
    std::vector<DataPointType> points_;

    // bounding box: the smallest rectangle containing all data DataPointType's. 
    RectangleType bbox_;

    // root node of the Kdtree. 
    Node *root_;

    // maximum number of points per leaf. 
    int leaf_nmax_;

  private:

    // helper functions
    void initialize();
    void merge_duplicates();
    RectangleType compute_bounding_box() const;
    Node* construct_tree(int, int, int);

    Node* tree_deep_copy(const Node*) const;

    void delete_tree(Node*);

    void retrieve_point_indices(const Node*, std::vector<IndexType>&) const;
    void retrieve_range_indices(const Node*, int, 
                                const RectangleType&, const RectangleType&, 
                                std::vector<IndexType>&) const;
    void retrieve_partitions(const Node *, int, int, 
                             const RectangleType&, 
                             std::vector<RectangleType>&) const;
};

// Implementations
// ---------------

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::print_partitions(int depth, std::ostream &os) const {
  std::vector<RectangleType> result;
  retrieve_partitions(root_, 0, depth, bbox_, result);
  for (const auto &r : result) { os << r << std::endl; }
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
inline void Kdtree<D,PAttrT,NAttrT,FloatT>::report_partitions(int depth, std::vector<RectangleType> &result) const {
  retrieve_partitions(root_, 0, depth, bbox_, result);
}

// use DFS to retrieve partitions. 
template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::retrieve_partitions(
    const Node *r, int d, int depth, 
    const RectangleType &partition,
    std::vector<RectangleType> &result) const {

  // handle null tree separately
  if (r == nullptr) { return; }

  // both conditions required. see API. 
  if (depth == 0 || r->is_leaf()) { result.push_back(partition); return; }

  retrieve_partitions(r->left_, (d+1)%D, depth-1, 
                      partition.lower_halfspace(d, r->split_), result);
  retrieve_partitions(r->right_, (d+1)%D, depth-1, 
                      partition.upper_halfspace(d, r->split_), result);

}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::print_range_search(const RectangleType &query_range, std::ostream &os) const {

  std::vector<IndexType> result_indices;
  retrieve_range_indices(root_, 0, bbox_, query_range, result_indices);
  for (auto i : result_indices) { os << points_[i] << std::endl; }

}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::range_search(const RectangleType &query_range, std::vector<DataPointType> &result) const {
  std::vector<IndexType> result_indices;
  retrieve_range_indices(root_, 0, bbox_, query_range, result_indices);
  for (auto i : result_indices) { result.emplace_back(points_[i]); }
}

// standard range search algorithm for kdtrees. 
template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::retrieve_range_indices(
    const Node *v, int d,
    const RectangleType &data_range, 
    const RectangleType &query_range, 
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
    RectangleType left_halfspace = data_range.lower_halfspace(d, v->split_);
    if (query_range.contains(left_halfspace)) {
      retrieve_point_indices(v->left_, result);
    } else {
      retrieve_range_indices(v->left_, (d+1)%D, left_halfspace, query_range, result);
    }

    // right halfspace
    RectangleType right_halfspace = data_range.upper_halfspace(d, v->split_);
    if (query_range.contains(right_halfspace)) {
      retrieve_point_indices(v->right_, result);
    } else {
      retrieve_range_indices(v->right_, (d+1)%D, right_halfspace, query_range, result);
    }

  }
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void swap(Kdtree<D,PAttrT,NAttrT,FloatT> &lhs, Kdtree<D,PAttrT,NAttrT,FloatT> &rhs) {

  using std::swap;
  
  swap(lhs.points_, rhs.points_);
  swap(lhs.bbox_, rhs.bbox_);
  swap(lhs.leaf_nmax_, rhs.leaf_nmax_);

  // simply switch the root_ pointers.
  typename Kdtree<D,PAttrT,NAttrT,FloatT>::Node *ptemp = lhs.root_;
  lhs.root_ = rhs.root_;
  rhs.root_ = ptemp;
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
inline Kdtree<D,PAttrT,NAttrT,FloatT>& Kdtree<D,PAttrT,NAttrT,FloatT>::operator=(Kdtree<D,PAttrT,NAttrT,FloatT> rhs) {
  swap(*this, rhs); return *this;
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
inline bool Kdtree<D,PAttrT,NAttrT,FloatT>::empty() const { return root_ == nullptr; }

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
inline typename Kdtree<D,PAttrT,NAttrT,FloatT>::IndexType Kdtree<D,PAttrT,NAttrT,FloatT>::size() const { return points_.size(); }

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
inline const typename Kdtree<D,PAttrT,NAttrT,FloatT>::RectangleType& 
Kdtree<D,PAttrT,NAttrT,FloatT>::bounding_box() const { return bbox_; }

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
inline int Kdtree<D,PAttrT,NAttrT,FloatT>::leaf_nmax() const 
{ return leaf_nmax_; }

// DFS to the leaves and print the index range of points to os
template<int D, typename PAttrT, typename NAttrT, typename FloatT>
inline void Kdtree<D,PAttrT,NAttrT,FloatT>::report_leaves(
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

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
inline const std::vector<typename Kdtree<D,PAttrT,NAttrT,FloatT>::DataPointType>&
Kdtree<D,PAttrT,NAttrT,FloatT>::points() const {
  return points_;
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::print_points(std::ostream &os) const {
  for (const auto &p : points_) { os << p << std::endl; }
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::retrieve_point_indices(const Node *r, std::vector<IndexType> &result) const {

  // check for null tree
  if (r == nullptr) { return; }

  // simply append all points under this node
  for (IndexType i = r->start_idx_; i <= r->end_idx_; ++i) {
    result.emplace_back(i);
  }
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
Kdtree<D,PAttrT,NAttrT,FloatT>::Kdtree() : points_(0), bbox_(), root_(nullptr) {}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
Kdtree<D,PAttrT,NAttrT,FloatT>::~Kdtree() { delete_tree(root_); }

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::initialize() {

  // merge duplicate keys
  merge_duplicates();

  // compute a minimum axis-aligned rectangle containing all the points
  bbox_ = compute_bounding_box();

  // build the tree
  if (!points_.empty()) { root_ = construct_tree(0, points_.size()-1, 0); }
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
Kdtree<D,PAttrT,NAttrT,FloatT>::Kdtree(
    const std::vector<DataPointType> &points, int leaf_nmax) 
  : points_(points), bbox_(), root_(nullptr), leaf_nmax_(leaf_nmax) {
  initialize();
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
Kdtree<D,PAttrT,NAttrT,FloatT>::Kdtree(
    std::vector<DataPointType> &&points, int leaf_nmax) 
  : points_(std::move(points)), bbox_(), root_(nullptr), leaf_nmax_(leaf_nmax) {
  initialize();
}

// merge duplicate keys in the data by merging the attributes. 
template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::merge_duplicates() {

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
template<int D, typename PAttrT, typename NAttrT, typename FloatT>
typename Kdtree<D,PAttrT,NAttrT,FloatT>::RectangleType Kdtree<D,PAttrT,NAttrT,FloatT>::compute_bounding_box() const {

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
// + this procedure rearranges points_ such that the indices stored at the leaves
//   refer to the appropriate data point. 
// Note: do NOT change the ordering of points_ outside of this function. 
template<int D, typename PAttrT, typename NAttrT, typename FloatT>
typename Kdtree<D,PAttrT,NAttrT,FloatT>::Node* Kdtree<D,PAttrT,NAttrT,FloatT>::construct_tree(int i, int j, int d) {

  Node *p = new Node();
  p->start_idx_ = i;
  p->end_idx_ = j;

  // create a leaf node when [i,j] contains no more than leaf_nmax_ indices. 
  if (j-i+1 <= leaf_nmax_) { 

    p->attr_ = NodeAttributesType::extract_point_attributes(points_[i]);
    for (IndexType k = i+1; k <= j; ++k) {
      p->attr_.merge(NodeAttributesType::extract_point_attributes(points_[k]));
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

    // pre-order tree walk. 
    p->left_ = construct_tree(i, m, (d+1)%D);
    p->right_ = construct_tree(m+1, j, (d+1)%D);

    p->attr_ = merge(p->left_->attr_, p->right_->attr_);
  }
  return p;
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
Kdtree<D,PAttrT,NAttrT,FloatT>::Kdtree(Kdtree<D,PAttrT,NAttrT,FloatT> &&obj) noexcept : 
  points_(std::move(obj.points_)), 
  bbox_(std::move(obj.bbox_)), 
  leaf_nmax_(std::move(obj.leaf_nmax_)) {
  root_ = obj.root_;
  obj.root_ = nullptr;
}

template<int D, typename PAttrT, typename NAttrT, typename FloatT>
Kdtree<D,PAttrT,NAttrT,FloatT>::Kdtree(const Kdtree<D,PAttrT,NAttrT,FloatT> &obj) 
  : points_(obj.points_), bbox_(obj.bbox_), leaf_nmax_(obj.leaf_nmax_) {
  root_ = tree_deep_copy(obj.root_);
}

// return a pointer to a deep copy of the tree pointed to by the argument. 
template<int D, typename PAttrT, typename NAttrT, typename FloatT>
typename Kdtree<D,PAttrT,NAttrT,FloatT>::Node* 
Kdtree<D,PAttrT,NAttrT,FloatT>::tree_deep_copy(
    const Kdtree<D,PAttrT,NAttrT,FloatT>::Node *target_r) const {
  if (target_r == nullptr) { return nullptr; }
  Node *r = new Node(*target_r);
  r->left_ = tree_deep_copy(target_r->left_);
  r->right_ = tree_deep_copy(target_r->right_);
  return r;
}

// delete Kdtree<>::Nodes by a post order tree walk. 
template<int D, typename PAttrT, typename NAttrT, typename FloatT>
void Kdtree<D,PAttrT,NAttrT,FloatT>::delete_tree(Node *p) {
  if (p == nullptr) { return; }
  if (!p->is_leaf()) { 
    delete_tree(p->left_);
    delete_tree(p->right_);
  }
  delete p; 
}

}

#endif
