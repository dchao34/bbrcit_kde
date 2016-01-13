#ifndef BBRCITKDE_KDTREE_H__
#define BBRCITKDE_KDTREE_H__

#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include <queue>
#include <iostream>

#include <Point.h>
#include <Rectangle.h>

// API
// ---

// TODO list: 
// + Handle duplicate keys. 

namespace bbrcit {

template <int D, typename T>
class Kdtree {

  using PointT = Point<D,T>;
  using IndexT = typename std::vector<PointT>::size_type;
  using RectangleT = Rectangle<D,T>;

  public: 
    Kdtree();
    Kdtree(const std::vector<PointT>&);
    ~Kdtree();

    bool empty() const;
    void report_leaves(std::ostream&) const;
    void report_leaves(std::vector<PointT>&) const;
    void range_search(std::vector<PointT>&, const RectangleT&) const;

  private:

    struct Node {
      T split_ = T();
      Node *left = nullptr, *right = nullptr;
      IndexT point_idx_ = 0;
      bool is_leaf() const { return left == nullptr && right == nullptr; }
    };

    std::vector<PointT> points_;
    RectangleT bbox_;
    Node *root_;

    Node* construct_tree(int, int, int);
    void delete_tree(Node*);

    RectangleT compute_bounding_box() const;
    void retrieve_leaves(const Node*, std::vector<IndexT>&) const;
    void retrieve_range(const Node*, int, 
                        const RectangleT&, const RectangleT&, 
                        std::vector<IndexT>&) const;
};

// Implementations
// ---------------

template<int D, typename T>
void Kdtree<D,T>::range_search(std::vector<PointT> &result, const RectangleT &query_range) const {

  std::vector<IndexT> result_indices;
  retrieve_range(root_, 0, bbox_, query_range, result_indices);

  result.reserve(result_indices.size());
  for (auto i : result_indices) { result.emplace_back(points_[i]); }
}

template<int D, typename T>
typename Kdtree<D,T>::Node* Kdtree<D,T>::construct_tree(int i, int j, int d) {
  Node *p = new Node();
  if (i == j) { 
    p->point_idx_ = i;
  } else {
    int m = i + (j-i) / 2;
    std::nth_element(points_.begin()+i, points_.begin()+m, points_.begin()+j+1, 
                    [d] (const PointT &p1, const PointT &p2) { return p1[d] < p2[d]; });
    p->split_ = points_[m][d];
    p->left = construct_tree(i, m, (d+1)%D);
    p->right = construct_tree(m+1, j, (d+1)%D);
  }
  return p;
}

template<int D, typename T>
void Kdtree<D,T>::retrieve_range(
    const Node *v, int d,
    const RectangleT &data_range, 
    const RectangleT &query_range, 
    std::vector<IndexT> &result) const {

  if (v == nullptr) { return; }

  if (v->is_leaf()) { 
    if (query_range.contains(points_[v->point_idx_])) {
      result.push_back(v->point_idx_);
    }
  } else {

    // left halfspace
    RectangleT left_halfspace = data_range.lower_halfspace(d, v->split_);
    if (query_range.contains(left_halfspace)) {
      retrieve_leaves(v->left, result);
    } else {
      retrieve_range(v->left, (d+1)%D, left_halfspace, query_range, result);
    }

    // right halfspace
    RectangleT right_halfspace = data_range.upper_halfspace(d, v->split_);
    if (query_range.contains(right_halfspace)) {
      retrieve_leaves(v->right, result);
    } else {
      retrieve_range(v->right, (d+1)%D, right_halfspace, query_range, result);
    }

  }
}

template<int D, typename T>
typename Kdtree<D,T>::RectangleT Kdtree<D,T>::compute_bounding_box() const {

  if (points_.empty()) { return RectangleT(); }

  std::vector<T> min_coord_val(D, std::numeric_limits<T>::max());
  std::vector<T> max_coord_val(D, std::numeric_limits<T>::min());
  for (const auto &p : points_) {
    for (int i = 0; i < D; ++i) { 
      min_coord_val[i] = std::min(p[i], min_coord_val[i]);
      max_coord_val[i] = std::max(p[i], max_coord_val[i]);
    }
  }

  RectangleT r;
  for (int i = 0; i < D; ++i) { 
    r.resize(i, {min_coord_val[i], max_coord_val[i]});
  }

  return r;
}

template<int D, typename T>
inline bool Kdtree<D,T>::empty() const { return root_ == nullptr; }

template<int D, typename T>
void Kdtree<D,T>::report_leaves(std::vector<PointT> &result) const {

  std::vector<IndexT> leaf_indices;
  retrieve_leaves(root_, leaf_indices);
  result.reserve(leaf_indices.size());
  for (auto i : leaf_indices) { result.emplace_back(points_[i]); }

}

template<int D, typename T>
void Kdtree<D,T>::report_leaves(std::ostream &os) const {

  std::vector<IndexT> leaf_indices;
  retrieve_leaves(root_, leaf_indices);
  for (auto i : leaf_indices) { os << points_[i] << std::endl; }

}

template<int D, typename T>
void Kdtree<D,T>::retrieve_leaves(const Node *r, std::vector<IndexT> &result) const {

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
      if (p->left) { q.push(p->left); }
      if (p->right) { q.push(p->right); }
    }
  }
}

template<int D, typename T>
Kdtree<D,T>::Kdtree() : points_(0), bbox_(), root_(nullptr) {}

template<int D, typename T>
Kdtree<D,T>::~Kdtree() { delete_tree(root_); }

template<int D, typename T>
Kdtree<D,T>::Kdtree(const std::vector<PointT> &points) : points_(points), bbox_(), root_(nullptr) {

  bbox_ = compute_bounding_box();

  if (!points_.empty()) {
    root_ = construct_tree(0, points_.size()-1, 0);
  }

}

template<int D, typename T>
void Kdtree<D,T>::delete_tree(Node *p) {
  if (p == nullptr) { return; }
  if (!p->is_leaf()) { 
    delete_tree(p->left);
    delete_tree(p->right);
  }
  delete p; 
}

}

#endif
