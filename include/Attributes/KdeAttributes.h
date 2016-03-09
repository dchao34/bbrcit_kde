#ifndef BBRCITKDE_KDEATTRIBUTES_H__
#define BBRCITKDE_KDEATTRIBUTES_H__

#include <limits>
#include <iostream>

#include <KdeTraits.h>

#include "AttributesInterface.h"

namespace bbrcit {

// KdeAttributes<>
// ---------------------

template<typename T> class KdeAttributes;

template<typename T> 
std::ostream& operator<<(std::ostream&, const KdeAttributes<T>&);

template<typename T>
class KdeAttributes {

  public: 

    KdeAttributes();
    KdeAttributes(T weight);

    KdeAttributes(const KdeAttributes<T>&) = default;
    KdeAttributes(KdeAttributes<T>&&) = default;
    KdeAttributes<T>& operator=(const KdeAttributes<T>&) = default;
    KdeAttributes<T>& operator=(KdeAttributes<T>&&) = default;
    ~KdeAttributes() = default;

    const T& weight() const { return weight_; }
    const T& lower() const { return lower_; }
    const T& upper() const { return upper_; }
    T middle() const { return lower_+(upper_-lower_)/2; }

    void set_weight(const T &w) { weight_ = w; }
    void set_lower(const T &l) { lower_ = l; }
    void set_upper(const T &u) { upper_ = u; }

    // required functions for the Attributes<> interface
    KdeAttributes<T>& merge(const KdeAttributes<T>&);

  private: 
    T weight_;
    T lower_, upper_;
};

template<typename T>
KdeAttributes<T>& 
KdeAttributes<T>::merge(const KdeAttributes<T> &rhs) {
  weight_ += rhs.weight_;
  lower_ = std::min(lower_, rhs.lower_);
  upper_ = std::max(upper_, rhs.upper_);
  return *this;
}

template<typename T> 
std::ostream& operator<<(std::ostream &os, const KdeAttributes<T> &att) {
  os << "("; 
  os << att.weight();
  os << ")";
  return os;
}

template<typename T>
KdeAttributes<T>::KdeAttributes() 
  : weight_(ConstantTraits<T>::one()),
    lower_(ConstantTraits<T>::zero()),
    upper_(ConstantTraits<T>::zero()) {
}

template<typename T>
KdeAttributes<T>::KdeAttributes(T w) 
  : weight_(w), 
    lower_(ConstantTraits<T>::zero()),
    upper_(ConstantTraits<T>::zero()) {
}

}

#endif
