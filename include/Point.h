#ifndef BBRCITKDE_POINT_H__
#define BBRCITKDE_POINT_H__

#include <vector>
#include <utility>
#include <exception>
#include <initializer_list>
#include <iostream>

// API
// ---

namespace bbrcit {

template <int D, typename T> 
class Point;

// custom swap for Point<>
template <int D, typename T> 
void swap(Point<D,T>&, Point<D,T>&);

// reads the next D space separated tokens into Point<D,T>.
template <int D, typename T> 
std::istream& operator>>(std::istream&, Point<D,T>&);

// prints Point<> as (coord1, coord2, ..., coordD).
template <int D, typename T> 
std::ostream& operator<<(std::ostream&, const Point<D,T>&);

// addition in euclidean space
template <int D, typename T> 
Point<D,T> operator+(const Point<D,T>&, const Point<D,T>&);

// subtraction in euclidean space
template <int D, typename T> 
Point<D,T> operator-(const Point<D,T>&, const Point<D,T>&);

// multiplication by scalar in euclidean space
template <int D, typename T> Point<D,T> operator*(const Point<D,T>&, double);
template <int D, typename T> Point<D,T> operator*(double, const Point<D,T>&);

// division by scalar in euclidean space
template <int D, typename T> Point<D,T> operator/(const Point<D,T>&, double);

// Point<> models a D dimensional point in euclidean space. 
template <int D, typename T = double>
class Point {

  friend void swap<>(Point<D,T>&, Point<D,T>&);

  public: 

    // default constructor: gives the zero element 
    Point();

    // constructor: allows the usage Point<2,double> p = { 0.0, 1.0 };
    Point(const std::initializer_list<T>&);

    // copy-control
    Point(const Point<D,T>&);
    Point<D,T>& operator=(const Point<D,T>&);
    Point(Point<D,T>&&) noexcept;
    Point<D,T>& operator=(Point<D,T> &&rhs) noexcept;
    ~Point();

    // reset coordinate values: coordinate i set to element i in arg. 
    // other elements set to T(). 
    void reset(const std::initializer_list<T> &li);

    // access coordinate values by indexing (0 for 1st coordiate etc.)
    const T& operator[](int) const;
    T& operator[](int);

    // the usual scaler multiplication and division in euclidean space
    Point<D,T>& operator*=(double);
    Point<D,T>& operator/=(double);

    // the usual addition and subtraction in euclidean space
    Point<D,T>& operator+=(const Point<D,T>&);
    Point<D,T>& operator-=(const Point<D,T>&);

  private:

    // represent the point using a list of coordinates of length D. 
    std::vector<T> coord_;
};


// Implementations
// ---------------

template <int D, typename T> 
void Point<D,T>::reset(const std::initializer_list<T> &li) {
  coord_ = li; 
  if (coord_.size() != D) { coord_.resize(D); }
}

template <int D, typename T> 
Point<D,T> operator/(const Point<D,T> &lhs, double rhs) {
  Point<D,T> result(lhs);
  result /= rhs;
  return result;
}

template <int D, typename T> 
Point<D,T> operator*(double lhs, const Point<D,T> &rhs) {
  Point<D,T> result(rhs);
  result *= lhs;
  return result;
}

template <int D, typename T> 
Point<D,T> operator*(const Point<D,T> &lhs, double rhs) {
  Point<D,T> result(lhs);
  result *= rhs;
  return result;
}

template <int D, typename T> 
Point<D,T> operator+(const Point<D,T> &lhs, const Point<D,T> &rhs) {
  Point<D,T> result(lhs);
  result += rhs;
  return result;
}

template <int D, typename T> 
Point<D,T> operator-(const Point<D,T> &lhs, const Point<D,T> &rhs) {
  Point<D,T> result(lhs);
  result -= rhs;
  return result;
}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator+=(const Point<D,T> &rhs) {
  for (int i = 0; i < D; ++i) { coord_[i] += (rhs.coord_)[i]; }
  return *this;
}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator-=(const Point<D,T> &rhs) {
  for (int i = 0; i < D; ++i) { coord_[i] -= (rhs.coord_)[i]; }
  return *this;
}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator*=(double c) {
  for (auto &e : coord_) { e *= c; }
  return *this;
}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator/=(double c) {
  if (!c) { 
    throw std::domain_error("Point<>: operator/=(): division by zero. ");
  }
  for (auto &e : coord_) { e /= c; }
  return *this;
}

template<int D, typename T>
inline void swap(Point<D,T> &lhs, Point<D,T> &rhs) {
  using std::swap; swap(lhs.coord_, rhs.coord_); return; 
}

template<int D, typename T>
Point<D,T>::Point() : coord_(D) {}

template<int D, typename T>
Point<D,T>::Point(const std::initializer_list<T> &li) : coord_(li) {
  if (coord_.size() != D) { coord_.resize(D); }
}

template<int D, typename T>
Point<D,T>::~Point() {}

template<int D, typename T>
Point<D,T>::Point(const Point<D,T> &rhs) : coord_(rhs.coord_) {}

template<int D, typename T>
Point<D,T>::Point(Point<D,T> &&t) noexcept : coord_(std::move(t.coord_)) {}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator=(const Point<D,T> &rhs) {
  if (this == &rhs) { return *this; }
  coord_ = rhs.coord_;
  return *this;
}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator=(Point<D,T> &&rhs) noexcept {
  if (this == &rhs) { return *this; }
  coord_ = std::move(rhs.coord_);
  return *this;
}

template<int D, typename T>
std::istream& operator>>(std::istream &is, Point<D, T> &p) {
  for (int i = 0; i < D && is; ++i) { is >> p[i]; }
  if (!is) { p = Point<D, T>(); }
  return is;
}

template<int D, typename T>
std::ostream& operator<<(std::ostream &os, const Point<D, T> &p) {
  os << "(";
  os << p[0]; for (int i = 1; i < D; ++i) { os << ", " << p[i]; }
  os << ")";
  return os;
}

template <int D, typename T>
inline const T& Point<D,T>::operator[](int i) const {
  return coord_.at(i);
}

template <int D, typename T>
inline T& Point<D,T>::operator[](int i) {
  return const_cast<T&>(static_cast<const Point<D,T>&>(*this)[i]);
}

}

#endif
