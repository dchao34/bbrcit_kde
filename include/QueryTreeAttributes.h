#ifndef BBRCITKDE_QUERYTREEATTRIBUTES_H__
#define BBRCITKDE_QUERYTREEATTRIBUTES_H__

#include <limits>
#include <iostream>

#include <AttributesInterface.h>
#include <KdeTraits.h>

namespace bbrcit {

// QueryTreeAttributes<>
// ---------------------

template<typename T> class QueryTreeAttributes;

template<typename T> 
std::ostream& operator<<(std::ostream&, const QueryTreeAttributes<T>&);

template<typename T>
class QueryTreeAttributes {

  public: 

    using WeightType = T;
    using BoundType = T;

    QueryTreeAttributes();

    QueryTreeAttributes(const T&, const T&, const T&);

    QueryTreeAttributes(const QueryTreeAttributes<T>&) = default;
    QueryTreeAttributes(QueryTreeAttributes<T>&&) = default;
    QueryTreeAttributes<T>& operator=(const QueryTreeAttributes<T>&) = default;
    QueryTreeAttributes<T>& operator=(QueryTreeAttributes<T>&&) = default;
    ~QueryTreeAttributes();

    const T& weight() const;
    void set_weight(const T&);

    const T& lower() const;
    T middle() const;
    const T& upper() const;
    void set_lower(const T&);
    void set_upper(const T&);

    // required functions for the Attributes<> interface
    QueryTreeAttributes<T>& merge(const QueryTreeAttributes<T>&);
    template<typename PointT> 
      static QueryTreeAttributes<T> extract_point_attributes(const PointT &p);

  private: 
    T weight_;
    T lower_, upper_;
};

template<typename T> 
  template<typename PointT> 
inline QueryTreeAttributes<T> 
QueryTreeAttributes<T>::extract_point_attributes(const PointT &p) {
  return p.attributes();
}

template<typename T>
QueryTreeAttributes<T>& 
QueryTreeAttributes<T>::merge(const QueryTreeAttributes<T> &rhs) {
  weight_ += rhs.weight_;
  lower_ = std::min(lower_, rhs.lower_);
  upper_ = std::max(upper_, rhs.upper_);
  return *this;
}

template<typename T> 
std::ostream& operator<<(std::ostream &os, const QueryTreeAttributes<T> &att) {
  os << "("; 
  os << att.weight() << ", ";
  os << att.lower() << ", ";
  os << att.upper();
  os << ")";
  return os;
}

template<typename T>
QueryTreeAttributes<T>::QueryTreeAttributes() 
  : weight_(ConstantTraits<T>::one()),
    lower_(ConstantTraits<T>::zero()),
    upper_(std::numeric_limits<T>::max()) {
}

template<typename T>
QueryTreeAttributes<T>::~QueryTreeAttributes() {}

template<typename T>
QueryTreeAttributes<T>::QueryTreeAttributes(const T &w, const T &l, const T &u) 
  : weight_(w), lower_(l), upper_(u) {
}

template<typename T>
inline const T& QueryTreeAttributes<T>::weight() const {
  return weight_;
}

template<typename T>
inline void QueryTreeAttributes<T>::set_weight(const T &w) {
  weight_ = w;
}

template<typename T>
inline const T& QueryTreeAttributes<T>::lower() const {
  return lower_;
}

template<typename T>
inline void QueryTreeAttributes<T>::set_lower(const T &l) {
  lower_ = l;
}

template<typename T>
inline const T& QueryTreeAttributes<T>::upper() const {
  return upper_;
}

template<typename T>
inline void QueryTreeAttributes<T>::set_upper(const T &u) {
  upper_ = u;
}

template<typename T>
inline T QueryTreeAttributes<T>::middle() const {
  return lower_ + (upper_ - lower_) / 2;
}

}

#endif
