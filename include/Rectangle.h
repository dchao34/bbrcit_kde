#ifndef BBRCITKDE_RECTANGLE_H__
#define BBRCITKDE_RECTANGLE_H__

#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>

#include <Point.h>
#include <Interval.h>

// API
// ---

namespace bbrcit {

template <int D, typename T> class Rectangle;

template <int D, typename T> 
void swap(Rectangle<D,T>&, Rectangle<D,T>&);

template <int D, typename T> 
std::ostream& operator<<(std::ostream&, const Rectangle<D,T>&);

template <int D, typename T> 
bool intersect(const Rectangle<D,T>&, const Rectangle<D,T>&);

template <int D, typename T=double>
class Rectangle {

  public: 

    using EdgeT = Interval<T>;
    friend void swap<>(Rectangle<D,T>&, Rectangle<D,T>&);

    Rectangle();
    Rectangle(const Point<D,T>&, const Point<D,T>&);
    Rectangle(const Rectangle<D,T>&);
    Rectangle(Rectangle<D,T>&&) noexcept;
    Rectangle& operator=(const Rectangle<D,T>&);
    Rectangle& operator=(Rectangle<D,T>&&) noexcept;
    ~Rectangle();

    // index access to edge intervals. 
    const EdgeT& operator[](int) const;
    EdgeT& operator[](int);

    // resize the interval of the d'th dimension to e.
    void resize(size_t d, const EdgeT &e);

    // resize this Rectangle<> to that defined by the Point<>'s in the argument.
    void resize(const Point<D,T>&, const Point<D,T>&);

    // returns true if the argument is fully contained in this Rectangle<>. 
    // specific choices for M include Point<>, Rectangle<>
    template <template <int D1, typename T1> class M> 
      bool contains(const M<D,T>&) const;

    // returns the (min, max) L2 distance from the Point<> to this Rectangle<>. 
    T min_dist(const Point<D,T>&) const;
    T max_dist(const Point<D,T>&) const;

  private:

    // Rectangle<D,T>'s are represented as a list of pair<T,T>'s of length D.
    // each pair<> represents the closed interval [lower, upper], where upper >= lower. 
    std::vector<EdgeT> intervals_;
};

// Implementations
// ---------------

template <int D, typename T>
void swap(Rectangle<D,T> &lhs, Rectangle<D,T> &rhs) {
  using std::swap; swap(lhs.intervals_, rhs.intervals_);
}

template <int D, typename T>
Rectangle<D,T>::Rectangle() : intervals_(D, EdgeT(T(), T()) ) {}

template <int D, typename T>
Rectangle<D,T>::~Rectangle() {}

template <int D, typename T>
Rectangle<D,T>::Rectangle(const Point<D,T> &p1, const Point<D,T> &p2) {
  for (int i = 0; i < D; ++i) { 
    auto p = std::minmax(p1[i], p2[i]);
    intervals_.emplace_back(p.first, p.second); 
  }
}

template <int D, typename T>
Rectangle<D,T>::Rectangle(const Rectangle &r) : intervals_(r.intervals_) {}

template <int D, typename T>
Rectangle<D,T>::Rectangle(Rectangle &&r) noexcept : intervals_(std::move(r.intervals_)) {}

template <int D, typename T>
Rectangle<D,T>& Rectangle<D,T>::operator=(const Rectangle &r) {
  if (this == &r) { return *this; }
  intervals_ = r.intervals_;
  return *this;
}

template <int D, typename T>
Rectangle<D,T>& Rectangle<D,T>::operator=(Rectangle &&r) noexcept {
  if (this == &r) { return *this; }
  intervals_ = std::move(r.intervals_);
  return *this;
}

template <int D, typename T>
T Rectangle<D,T>::min_dist(const Point<D,T> &p) const {
  T total = T(), curr = T();
  for (int i = 0; i < D; ++i) {
    curr = intervals_[i].min_dist(p[i]);
    curr *= curr;
    total += curr;
  }
  return std::sqrt(total);
}

template <int D, typename T>
T Rectangle<D,T>::max_dist(const Point<D,T> &p) const {
  T total = T(), curr = T();
  for (int i = 0; i < D; ++i) {
    curr = intervals_[i].max_dist(p[i]);
    curr *= curr;
    total += curr;
  }
  return std::sqrt(total);
}

template <int D, typename T>
  template <template <int D1, typename T1> class M>
bool Rectangle<D,T>::contains(const M<D,T> &obj) const {
  bool is_contained = true;
  for (int i = 0; i < D && is_contained; ++i) {
    is_contained = intervals_[i].contains(obj[i]);
  }
  return is_contained;
}

template <int D, typename T>
inline void Rectangle<D,T>::resize(size_t d, const EdgeT &e) { (*this)[d] = e; }

template <int D, typename T>
void Rectangle<D,T>::resize(const Point<D,T> &p1, const Point<D,T> &p2) { 
  for (int i = 0; i < D; ++i) { 
    auto p = std::minmax(p1[i], p2[i]);
    intervals_[i] = EdgeT(p.first, p.second);
  }
}

template <int D, typename T> 
bool intersect(const Rectangle<D,T> &lhs, const Rectangle<D,T> &rhs) {
  bool is_intersect = true;
  for (int i = 0; i < D && is_intersect; ++i) { 
    is_intersect = intersect(lhs[i], rhs[i]); 
  }
  return is_intersect;
}

template <int D, typename T> 
std::ostream& operator<<(std::ostream &os, const Rectangle<D,T> &r) {
  os << "{ "; os << r[0];
  for (int i = 1; i < D; ++i) { os << ", "; os << r[i]; }
  os << " }"; 
  return os;
}

template <int D, typename T>
inline const typename Rectangle<D,T>::EdgeT& 
Rectangle<D,T>::operator[](int i) const { return intervals_.at(i); }

template <int D, typename T>
inline typename Rectangle<D,T>::EdgeT& 
Rectangle<D,T>::operator[](int i) { 
  return const_cast<EdgeT&>(static_cast<const Rectangle&>(*this)[i]);
}


}

#endif
