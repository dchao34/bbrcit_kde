#ifndef BBRCITKDE_KDTREE_H__
#define BBRCITKDE_KDTREE_H__

#include <vector>
#include <memory>
#include <algorithm>
#include <queue>
#include <iostream>

#include <Point.h>

// API
// ---

namespace bbrcit {

template <int D, typename T>
class Kdtree {

  // TODO: Handle duplicate keys. 

  using point_type = Point<D,T>;

  public: 
    Kdtree();
    Kdtree(const std::vector<point_type>&);
    ~Kdtree();

    void print() const;

  private:

    struct Node {
      T split_ = T();
      Node *left = nullptr, *right = nullptr;
      size_t point_idx_ = 0;
      bool is_leaf() const { return left == nullptr && right == nullptr; }
    };

    std::vector<point_type> points_;
    Node *root_;

    Node* construct_tree(int, int, int);
    void delete_tree(Node*);
};

// Implementations
// ---------------

template<int D, typename T>
void Kdtree<D,T>::print() const {
  std::queue<const Node*> q;
  q.push(root_);
  while (!q.empty()) {
    auto p = q.front(); q.pop();
    if (p->is_leaf()) { continue; }
    std::cout << p->split_ << " ";
    if (p->left) { q.push(p->left); }
    if (p->right) { q.push(p->right); }
  }
}

template<int D, typename T>
Kdtree<D,T>::Kdtree() : points_(0), root_(nullptr) {}

template<int D, typename T>
Kdtree<D,T>::~Kdtree() { delete_tree(root_); }

template<int D, typename T>
Kdtree<D,T>::Kdtree(const std::vector<point_type> &points) : points_(points), root_(nullptr) {
  if (!points_.empty()) {
    root_ = construct_tree(0, points_.size()-1, 0);
  }
}

template<int D, typename T>
typename Kdtree<D,T>::Node* Kdtree<D,T>::construct_tree(int i, int j, int d) {
  Node *p = new Node();
  if (i == j) { 
    p->point_idx_ = i;
  } else {
    int m = i + (j-i) / 2;
    std::nth_element(points_.begin()+i, points_.begin()+m, points_.end(), 
                    [d] (const point_type &p1, const point_type &p2) { return p1[d] < p2[d]; });
    p->split_ = points_[m][d];
    p->left = construct_tree(i, m, (d+1)%D);
    p->right = construct_tree(m+1, j, (d+1)%D);
  }
  return p;
}

template<int D, typename T>
void Kdtree<D,T>::delete_tree(Node *p) {
  if (p == nullptr) { return; }
  if (p->is_leaf()) { 
    delete p; 
  } else {
    delete_tree(p->left);
    delete_tree(p->right);
  }
}

}

#endif
