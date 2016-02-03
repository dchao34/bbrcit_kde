#ifndef BBRCITKDE_POINTWEIGHTS_H__
#define BBRCITKDE_POINTWEIGHTS_H__

#include <iostream>

#include <KdeTraits.h>

#include "AttributesInterface.h"

namespace bbrcit {

// PointWeights<>
// ------------------

template<typename T> class PointWeights;

// prints PointWeights<> as { (weight) } to os. 
template<typename T> 
std::istream& operator>>(std::istream&, PointWeights<T>&);

// reads in a single token
template<typename T> 
std::ostream& operator<<(std::ostream&, const PointWeights<T>&);

// PointWeights<> contains only a weight attribute for Kdtree data points. 
template<typename T>
class PointWeights {

  friend std::istream& operator>><>(std::istream&, PointWeights<T>&);

  public: 

    using WeightType = T;

    // default constructor: weight equal to 1.
    PointWeights();

    // one-argument constructor: user configured weight. 
    PointWeights(const T&);

    PointWeights(const PointWeights<T>&) = default;
    PointWeights(PointWeights<T>&&) = default;
    PointWeights<T>& operator=(const PointWeights<T>&) = default;
    PointWeights<T>& operator=(PointWeights<T>&&) = default;
    ~PointWeights();

    const T& weight() const;
    void set_weight(const T&);

    // required functions for the Attributes<> interface
    PointWeights<T>& merge(const PointWeights<T>&);

  private: 
    T weight_;
};

template<typename T>
PointWeights<T>& 
PointWeights<T>::merge(const PointWeights<T> &rhs) {
  weight_ += rhs.weight_;
  return *this;
}

template<typename T> 
std::ostream& operator<<(std::ostream &os, const PointWeights<T> &att) {
  os << "(" << att.weight() << ")";
  return os;
}

template<typename T>
std::istream& operator>>(std::istream &is, PointWeights<T> &att) {
  is >> att.weight_;
  if (!is) { att.weight_ = ConstantTraits<T>::one(); }
  return is;
}

template<typename T>
PointWeights<T>::PointWeights() 
  : weight_(ConstantTraits<T>::one()) {
}

template<typename T>
PointWeights<T>::~PointWeights() {}

template<typename T>
PointWeights<T>::PointWeights(const T &w) 
  : weight_(w) {
}

template<typename T>
inline const T& PointWeights<T>::weight() const {
  return weight_;
}

template<typename T>
inline void PointWeights<T>::set_weight(const T &w) {
  weight_ = w;
}

}

#endif
