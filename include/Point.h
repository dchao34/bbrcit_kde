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

template <int D, typename T> 
void swap(Point<D,T>&, Point<D,T>&);

template <int D, typename T> 
std::istream& operator>>(std::istream&, Point<D,T>&);

template <int D, typename T> 
std::ostream& operator>>(std::ostream&, const Point<D,T>&);

template <int D, typename T> 
Point<D,T> operator+(const Point<D,T>&, const Point<D,T>&);

template <int D, typename T> 
Point<D,T> operator-(const Point<D,T>&, const Point<D,T>&);

template <int D, typename T> 
Point<D,T> operator*(const Point<D,T>&, double);

template <int D, typename T> 
Point<D,T> operator*(double, const Point<D,T>&);

template <int D, typename T> 
Point<D,T> operator/(const Point<D,T>&, double);

template <int D, typename T = double>
class Point {

  friend void swap<>(Point<D,T>&, Point<D,T>&);

  public: 

    // constructors
    Point();
    Point(const std::initializer_list<T>&);

    // copy-control
    Point(const Point<D,T>&);
    Point<D,T>& operator=(const Point<D,T>&);
    Point(Point<D,T>&&) noexcept;
    Point<D,T>& operator=(Point<D,T> &&rhs) noexcept;
    ~Point();

    // access coordinate values by indexing. 
    const T& operator[](int) const;
    T& operator[](int);

    // scaler multiplication and division
    Point<D,T>& operator*=(double);
    Point<D,T>& operator/=(double);

    // vector addition and subtraction
    Point<D,T>& operator+=(const Point<D,T>&);
    Point<D,T>& operator-=(const Point<D,T>&);

  private:
    std::vector<T> coord;
};


// Implementations
// ---------------

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
  for (int i = 0; i < D; ++i) { coord[i] += (rhs.coord)[i]; }
  return *this;
}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator-=(const Point<D,T> &rhs) {
  for (int i = 0; i < D; ++i) { coord[i] -= (rhs.coord)[i]; }
  return *this;
}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator*=(double c) {
  for (auto &e : coord) { e *= c; }
  return *this;
}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator/=(double c) {
  if (!c) { 
    throw std::domain_error("Point<>: operator/=(): division by zero. ");
  }
  for (auto &e : coord) { e /= c; }
  return *this;
}

template<int D, typename T>
inline void swap(Point<D,T> &lhs, Point<D,T> &rhs) {
  using std::swap; swap(lhs.coord, rhs.coord); return; 
}

template<int D, typename T>
Point<D,T>::Point() : coord(D) {}

template<int D, typename T>
Point<D,T>::Point(const std::initializer_list<T> &li) : coord(li) {}

template<int D, typename T>
Point<D,T>::~Point() {}

template<int D, typename T>
Point<D,T>::Point(const Point<D,T> &rhs) : coord(rhs.coord) {}

template<int D, typename T>
Point<D,T>::Point(Point<D,T> &&t) noexcept : coord(std::move(t.coord)) {}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator=(const Point<D,T> &rhs) {
  if (this == &rhs) { return *this; }
  coord = rhs.coord;
  return *this;
}

template<int D, typename T>
Point<D,T>& Point<D,T>::operator=(Point<D,T> &&rhs) noexcept {
  if (this == &rhs) { return *this; }
  coord = std::move(rhs.coord);
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
  os << p[0]; for (int i = 1; i < D; ++i) { os << " " << p[i]; }
  return os;
}

template <int D, typename T>
inline const T& Point<D,T>::operator[](int i) const {
  return coord.at(i);
}

template <int D, typename T>
inline T& Point<D,T>::operator[](int i) {
  return const_cast<T&>(static_cast<const Point<D,T>&>(*this)[i]);
}

}

#endif
