#ifndef BBRCITKDE_KDTREE_H__
#define BBRCITKDE_KDTREE_H__

#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include <queue>
#include <iostream>

#include <DecoratedPoint.h>
#include <Rectangle.h>
#include <FloatUtils.h>
#include <PointWeights.h>

// API
// ---

namespace bbrcit {

template<int D, typename AttrT, typename FloatT> class Kdtree;

template<int D, typename AttrT, typename FloatT>
void swap(Kdtree<D,AttrT,FloatT>&, Kdtree<D,AttrT,FloatT>&);

// Kdtree<> implements a D-dimensional kdtree. It currently supports:
// + Range search. 
template<int D, 
         typename AttrT=PointWeights<int>, 
         typename FloatT = double>
class Kdtree {

  public: 

    using PointType = DecoratedPoint<D,AttrT,FloatT>;
    using FloatType = FloatT;
    using RectangleType = Rectangle<D,FloatType>;
    static int dim() { return D; }

  private:
    using IndexType = typename std::vector<PointType>::size_type;
    friend void swap<>(Kdtree<D,AttrT,FloatT>&, Kdtree<D,AttrT,FloatT>&);

  public: 

    // default constructor yields a null tree. 
    Kdtree();

    // construct Kdtree out of data given as a list of PointType's. 
    Kdtree(const std::vector<PointType> &data);
    Kdtree(std::vector<PointType> &&data);

    // copy-control. operator='s use copy and swap.
    Kdtree(const Kdtree<D,AttrT,FloatT>&);
    Kdtree(Kdtree<D,AttrT,FloatT>&&) noexcept;
    Kdtree<D,AttrT,FloatT>& operator=(Kdtree<D,AttrT,FloatT>);
    ~Kdtree();

    // returns true if this is a null tree. 
    bool empty() const;

    // returns the number of data points in this Kdtree.
    IndexType size() const;

    // returns the bounding box of the data points.
    RectangleType get_bounding_box() const;

    // (1) prints each PointType stored in the Kdtree on a separate line os.
    // (2) returns the PointType's stored in the Kdtree in a vector.
    void report_leaves(std::ostream &os) const;
    void report_leaves(std::vector<PointType>&) const;

    // (1) print each partition at depth d and leaf partitions of depth <= d on a separate line of os.
    // (2) returns each partition at depth d and leaf partitions of depth <= d in a vector.
    void report_partitions(int d, std::ostream&) const;
    void report_partitions(int d, std::vector<RectangleType>&) const;

    // (1) print each PointType contained in the query window on a separate line of os.
    // (2) returns each PointType contained in the query window in a vector. 
    void range_search(const RectangleType &query, std::ostream &os) const;
    void range_search(const RectangleType &query, std::vector<PointType>&) const;

  private:

    // Kdtree<>::Node represents a node in the Kdtree. 
    // (I/L) are members that are meaningful for internal/leaf nodes. 
    // + An object represents an internal node iff left_=right_=nullptr.
    // + Leaf nodes store an index (point_idx_) into the array of data PointType's. 
    struct Node {

      // (I) coordinate at which to partition half spaces.
      FloatType split_ = FloatType();

      // (I/L) links to the daughters
      Node *left_ = nullptr, *right_ = nullptr;

      // (L) links to the daughters
      IndexType point_idx_ = 0;

      // returns true if this object is a leaf
      bool is_leaf() const { return left_ == nullptr && right_ == nullptr; }

    };

    // list of PointType's that stores the data from which to construct the Kdtree. 
    // Note: the order of PointType's should NOT change once the object is constructed.
    std::vector<PointType> points_;

    // bounding box: the smallest rectangle containing all data PointType's. 
    RectangleType bbox_;

    // root node of the Kdtree. 
    Node *root_;

    // helper functions
    void initialize();
    void merge_duplicates();
    RectangleType compute_bounding_box() const;
    Node* construct_tree(int, int, int);

    Node* tree_deep_copy(const Node*) const;

    void delete_tree(Node*);

    void retrieve_leaves(const Node*, std::vector<IndexType>&) const;
    void retrieve_range(const Node*, int, 
                        const RectangleType&, const RectangleType&, 
                        std::vector<IndexType>&) const;
    void retrieve_partitions(const Node *, int, int, 
                             const RectangleType&, 
                             std::vector<RectangleType>&) const;
};

// Implementations
// ---------------

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::report_partitions(int depth, std::ostream &os) const {
  std::vector<RectangleType> result;
  retrieve_partitions(root_, 0, depth, bbox_, result);
  for (const auto &r : result) { os << r << std::endl; }
}

template<int D, typename AttrT, typename FloatT>
inline void Kdtree<D,AttrT,FloatT>::report_partitions(int depth, std::vector<RectangleType> &result) const {
  retrieve_partitions(root_, 0, depth, bbox_, result);
}

// use DFS to retrieve partitions. 
template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::retrieve_partitions(
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

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::range_search(const RectangleType &query_range, std::ostream &os) const {

  std::vector<IndexType> result_indices;
  retrieve_range(root_, 0, bbox_, query_range, result_indices);
  for (auto i : result_indices) { os << points_[i] << std::endl; }

}

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::range_search(const RectangleType &query_range, std::vector<PointType> &result) const {

  std::vector<IndexType> result_indices;
  retrieve_range(root_, 0, bbox_, query_range, result_indices);

  result.reserve(result_indices.size());
  for (auto i : result_indices) { result.emplace_back(points_[i]); }
}

// standard range search algorithm for kdtrees. 
template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::retrieve_range(
    const Node *v, int d,
    const RectangleType &data_range, 
    const RectangleType &query_range, 
    std::vector<IndexType> &result) const {

  if (v == nullptr) { return; }

  if (v->is_leaf()) { 
    if (query_range.contains(points_[v->point_idx_])) {
      result.push_back(v->point_idx_);
    }
  } else {

    // left halfspace
    RectangleType left_halfspace = data_range.lower_halfspace(d, v->split_);
    if (query_range.contains(left_halfspace)) {
      retrieve_leaves(v->left_, result);
    } else {
      retrieve_range(v->left_, (d+1)%D, left_halfspace, query_range, result);
    }

    // right halfspace
    RectangleType right_halfspace = data_range.upper_halfspace(d, v->split_);
    if (query_range.contains(right_halfspace)) {
      retrieve_leaves(v->right_, result);
    } else {
      retrieve_range(v->right_, (d+1)%D, right_halfspace, query_range, result);
    }

  }
}

template<int D, typename AttrT, typename FloatT>
void swap(Kdtree<D,AttrT,FloatT> &lhs, Kdtree<D,AttrT,FloatT> &rhs) {

  using std::swap;
  
  swap(lhs.points_, rhs.points_);
  swap(lhs.bbox_, rhs.bbox_);

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
inline typename Kdtree<D,AttrT,FloatT>::IndexType Kdtree<D,AttrT,FloatT>::size() const { return points_.size(); }

template<int D, typename AttrT, typename FloatT>
inline typename Kdtree<D,AttrT,FloatT>::RectangleType Kdtree<D,AttrT,FloatT>::get_bounding_box() const { return bbox_; }

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::report_leaves(std::vector<PointType> &result) const {

  std::vector<IndexType> leaf_indices;
  retrieve_leaves(root_, leaf_indices);
  result.reserve(leaf_indices.size());
  for (auto i : leaf_indices) { result.emplace_back(points_[i]); }

}

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::report_leaves(std::ostream &os) const {

  std::vector<IndexType> leaf_indices;
  retrieve_leaves(root_, leaf_indices);
  for (auto i : leaf_indices) { os << points_[i] << std::endl; }

}

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::retrieve_leaves(const Node *r, std::vector<IndexType> &result) const {

  // check for null tree
  if (r == nullptr) { return; }

  // BFS for the leaves
  std::queue<const Node*> q;
  q.push(r);
  while (!q.empty()) {
    auto p = q.front(); q.pop();
    if (p->is_leaf()) {
      result.emplace_back(p->point_idx_);
    } else {
      if (p->left_) { q.push(p->left_); }
      if (p->right_) { q.push(p->right_); }
    }
  }
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree() : points_(0), bbox_(), root_(nullptr) {}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::~Kdtree() { delete_tree(root_); }

template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::initialize() {

  // merge duplicate keys
  merge_duplicates();

  // compute a minimum axis-aligned rectangle containing all the points
  bbox_ = compute_bounding_box();

  // build the tree
  if (!points_.empty()) { root_ = construct_tree(0, points_.size()-1, 0); }
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree(const std::vector<PointType> &points) : points_(points), bbox_(), root_(nullptr) {
  initialize();
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree(std::vector<PointType> &&points) : points_(std::move(points)), bbox_(), root_(nullptr) {
  initialize();
}

// merge duplicate keys in the data by merging the attributes. 
template<int D, typename AttrT, typename FloatT>
void Kdtree<D,AttrT,FloatT>::merge_duplicates() {

  if (points_.empty()) { return; }

  // preprocess by sorting the data points lexicographically
  // note: exact comparison of floating point is ok for this purpose
  std::sort(points_.begin(), points_.end(), ExactLexicoLess<PointType>);

  // remove duplicates by merging attributes. 
  // algorithm is similar to the partition step in quicksort.
  size_t i = 0;
  for (size_t j = 1; j < points_.size(); ++j) {
    if (!ExactEqual(points_[i], points_[j])) {
      swap(points_[++i], points_[j]);
    } else {
      points_[i].set_attributes(
          merge(points_[i].get_attributes(), 
                points_[j].get_attributes())
      );
    }
  }

  // with duplicates removed, ok to shrink data size
  points_.resize(i+1);

}


// if points_ is non-empty, return the minimum area axis-aligned rectangle that 
// containing all PointType's in points_. otherwise return a degenerate rectangle 
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
// + this procedure rearranges points_ such that the indices stored at the leaves
//   refer to the appropriate data point. 
// Note: do NOT change the ordering of points_ outside of this function. 
template<int D, typename AttrT, typename FloatT>
typename Kdtree<D,AttrT,FloatT>::Node* Kdtree<D,AttrT,FloatT>::construct_tree(int i, int j, int d) {

  Node *p = new Node();

  // create a leaf node when there is exactly one data point in [i,j].
  // feature suggestion: allow n_max data points at leaves. 
  if (i == j) { 
    p->point_idx_ = i;

  // partition points_ by the lower median m. [i,m] go in the left subtree while 
  // [m+1,j] go in the right subtree. 
  } else {

    // median finding in expected linear time. Note: C++ requires the end iterator
    // to be one past the last index, hence the +1 in the 3rd argument. 
    int m = i + (j-i) / 2;
    std::nth_element(points_.begin()+i, points_.begin()+m, points_.begin()+j+1, 
                    [d] (const PointType &p1, const PointType &p2) { return p1[d] < p2[d]; });
    p->split_ = points_[m][d];

    // pre-order tree walk. 
    p->left_ = construct_tree(i, m, (d+1)%D);
    p->right_ = construct_tree(m+1, j, (d+1)%D);
  }
  return p;
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree(Kdtree<D,AttrT,FloatT> &&obj) noexcept : 
  points_(std::move(obj.points_)), 
  bbox_(std::move(obj.bbox_)) {
  root_ = obj.root_;
  obj.root_ = nullptr;
}

template<int D, typename AttrT, typename FloatT>
Kdtree<D,AttrT,FloatT>::Kdtree(const Kdtree<D,AttrT,FloatT> &obj) 
  : points_(obj.points_), bbox_(obj.bbox_) {
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
