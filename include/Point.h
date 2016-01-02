#ifndef BBRCITKDE_POINT_H__
#define BBRCITKDE_POINT_H__

#include <array>
#include <utility>
#include <iostream>

template <int D, typename T = double>
class Point {

  public: 

    Point() = default;
    template <typename ...E> Point(E&&... elems);

    ~Point() = default;
    Point(const Point<D,T>&) = default;

    Point<D,T>& operator=(const Point<D,T>&) = default;

    const T& operator[](int) const;
    T& operator[](int);

  private:
    std::array<T, D> coord;
};

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

// Is there a better solution for list initialization? 
// (current implementation based on stack overflow question 6893700)
// Forwarding is great, but brace initializers forbids narrowing; 
// e.g. can't do Point<3, double> = { 1, 2.0, 3.0 }.
template <int D, typename T>
template <typename ...E> 
Point<D,T>::Point(E&&... elems) : coord{{std::forward<E>(elems)...}} {};

template <int D, typename T>
inline const T& Point<D,T>::operator[](int i) const {
  return coord.at(i);
}

template <int D, typename T>
inline T& Point<D,T>::operator[](int i) {
  return const_cast<T&>(static_cast<const Point<D,T>&>(*this)[i]);
}

#endif
