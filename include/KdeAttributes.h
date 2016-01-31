#ifndef BBRCITKDE_KDEATTRIBUTES_H__
#define BBRCITKDE_KDEATTRIBUTES_H__

#include <limits>
#include <iostream>

#include <AttributesInterface.h>
#include <KdeTraits.h>

namespace bbrcit {

// KdeAttributes<>
// ---------------------

template<typename T> class KdeAttributes;

template<typename T> 
std::ostream& operator<<(std::ostream&, const KdeAttributes<T>&);

template<typename T>
class KdeAttributes {

  public: 

    using WeightType = T;
    using BoundType = T;

    KdeAttributes();

    KdeAttributes(const T&, const T&, const T&);

    KdeAttributes(const KdeAttributes<T>&) = default;
    KdeAttributes(KdeAttributes<T>&&) = default;
    KdeAttributes<T>& operator=(const KdeAttributes<T>&) = default;
    KdeAttributes<T>& operator=(KdeAttributes<T>&&) = default;
    ~KdeAttributes();

    const T& weight() const;
    void set_weight(const T&);

    const T& lower() const;
    T middle() const;
    const T& upper() const;
    void set_lower(const T&);
    void set_upper(const T&);

    // required functions for the Attributes<> interface
    KdeAttributes<T>& merge(const KdeAttributes<T>&);
    template<typename PointT> 
      static KdeAttributes<T> extract_point_attributes(const PointT &p);

  private: 
    T weight_;
    T lower_, upper_;
};

template<typename T> 
  template<typename PointT> 
inline KdeAttributes<T> 
KdeAttributes<T>::extract_point_attributes(const PointT &p) {
  return p.attributes();
}

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
  os << att.weight() << ", ";
  os << att.lower() << ", ";
  os << att.upper();
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
KdeAttributes<T>::~KdeAttributes() {}

template<typename T>
KdeAttributes<T>::KdeAttributes(const T &w, const T &l, const T &u) 
  : weight_(w), lower_(l), upper_(u) {
}

template<typename T>
inline const T& KdeAttributes<T>::weight() const {
  return weight_;
}

template<typename T>
inline void KdeAttributes<T>::set_weight(const T &w) {
  weight_ = w;
}

template<typename T>
inline const T& KdeAttributes<T>::lower() const {
  return lower_;
}

template<typename T>
inline void KdeAttributes<T>::set_lower(const T &l) {
  lower_ = l;
}

template<typename T>
inline const T& KdeAttributes<T>::upper() const {
  return upper_;
}

template<typename T>
inline void KdeAttributes<T>::set_upper(const T &u) {
  upper_ = u;
}

template<typename T>
inline T KdeAttributes<T>::middle() const {
  return lower_ + (upper_ - lower_) / 2;
}

}

#endif
