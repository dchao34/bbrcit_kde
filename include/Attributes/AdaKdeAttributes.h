#ifndef BBRCITKDE_ADAKDEATTRIBUTES_H__
#define BBRCITKDE_ADAKDEATTRIBUTES_H__

#include <limits>
#include <iostream>

#include <KdeTraits.h>

#include "AttributesInterface.h"

namespace bbrcit {

// AdaKdeAttributes<>
// ---------------------

template<typename T> class AdaKdeAttributes;

template<typename T> 
std::ostream& operator<<(std::ostream&, const AdaKdeAttributes<T>&);

template<typename T>
class AdaKdeAttributes {

  public: 

    AdaKdeAttributes();
    AdaKdeAttributes(T weight);

    AdaKdeAttributes(const AdaKdeAttributes<T>&) = default;
    AdaKdeAttributes(AdaKdeAttributes<T>&&) = default;
    AdaKdeAttributes<T>& operator=(const AdaKdeAttributes<T>&) = default;
    AdaKdeAttributes<T>& operator=(AdaKdeAttributes<T>&&) = default;
    ~AdaKdeAttributes() = default;

    // point weight
    const T& weight() const { return weight_; }
    void set_weight(const T &w) { weight_ = w; }

    // point mass
    const T& mass() const { return mass_; }
    void set_mass(const T &m) { mass_ = m; }

    // adaptive bandwidth correction and upper/lower bounds
    T abw() const { return lower_abw_+(upper_abw_-lower_abw_)/2; }
    const T& lower_abw() const { return lower_abw_; }
    const T& upper_abw() const { return upper_abw_; }
    void set_lower_abw(const T &l) { lower_abw_ = l; }
    void set_upper_abw(const T &u) { upper_abw_ = u; }

    // estimated density value and upper/lower bounds
    T value() const { return lower_+(upper_-lower_)/2; }
    const T& lower() const { return lower_; }
    const T& upper() const { return upper_; }
    void set_lower(const T &l) { lower_ = l; }
    void set_upper(const T &u) { upper_ = u; }


    // required functions for the Attributes<> interface
    AdaKdeAttributes<T>& merge(const AdaKdeAttributes<T>&);

  private: 
    T weight_;
    T mass_;
    T lower_abw_, upper_abw_;
    T lower_, upper_;
};

template<typename T>
AdaKdeAttributes<T>& 
AdaKdeAttributes<T>::merge(const AdaKdeAttributes<T> &rhs) {
  weight_ += rhs.weight_;
  mass_ += rhs.mass_;
  lower_abw_ = std::min(lower_abw_, rhs.lower_abw_);
  upper_abw_ = std::max(upper_abw_, rhs.upper_abw_);
  lower_ = std::min(lower_, rhs.lower_);
  upper_ = std::max(upper_, rhs.upper_);
  return *this;
}

template<typename T> 
std::ostream& operator<<(std::ostream &os, const AdaKdeAttributes<T> &att) {
  os << "("; 
  os << att.weight();
  os << ")";
  return os;
}

template<typename T>
AdaKdeAttributes<T>::AdaKdeAttributes() 
  : weight_(ConstantTraits<T>::one()),
    mass_(ConstantTraits<T>::one()),
    lower_abw_(ConstantTraits<T>::one()),
    upper_abw_(ConstantTraits<T>::one()),
    lower_(ConstantTraits<T>::zero()),
    upper_(ConstantTraits<T>::zero()) {
}

template<typename T>
AdaKdeAttributes<T>::AdaKdeAttributes(T w) 
  : weight_(w), 
    mass_(w),
    lower_abw_(ConstantTraits<T>::one()),
    upper_abw_(ConstantTraits<T>::one()),
    lower_(ConstantTraits<T>::zero()),
    upper_(ConstantTraits<T>::zero()) {
}

}

#endif
