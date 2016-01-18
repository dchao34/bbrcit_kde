#ifndef BBRCITKDE_POINTATTRIBUTES_H__
#define BBRCITKDE_POINTATTRIBUTES_H__

#include <iostream>
#include <KdeTraits.h>

namespace bbrcit {

// Attribute type definitions
// --------------------------

template<typename WeightT> class MinimalAttributes;

// prints MinimalAttributes<> as { (weight) } to os. 
template<typename WeightT> 
std::istream& operator>>(std::istream&, MinimalAttributes<WeightT>&);

// reads in a single token
template<typename WeightT> 
std::ostream& operator<<(std::ostream&, const MinimalAttributes<WeightT>&);

// MinimalAttributes<> contains the minimal attributes required of a Kdtree data point. 
template<typename WeightT>
class MinimalAttributes {

  friend std::istream& operator>><>(std::istream&, MinimalAttributes<WeightT>&);

  public: 

    using WeightType = WeightT;

    // default constructor: weight equal to 1.
    MinimalAttributes();

    // one-argument constructor: user configured weight. 
    MinimalAttributes(const WeightT&);

    MinimalAttributes(const MinimalAttributes<WeightT>&) = default;
    MinimalAttributes(MinimalAttributes<WeightT>&&) = default;
    MinimalAttributes<WeightT>& operator=(const MinimalAttributes<WeightT>&) = default;
    MinimalAttributes<WeightT>& operator=(MinimalAttributes<WeightT>&&) = default;
    ~MinimalAttributes();

    const WeightT& get_weight() const;
    void set_weight(const WeightT&);

    MinimalAttributes<WeightT>& merge(const MinimalAttributes<WeightT>&);

  private: 
    WeightT weight_;
};

template<typename WeightT>
MinimalAttributes<WeightT>& 
MinimalAttributes<WeightT>::merge(const MinimalAttributes<WeightT> &rhs) {
  weight_ += rhs.weight_;
  return *this;
}

template<typename WeightT> 
std::ostream& operator<<(std::ostream &os, const MinimalAttributes<WeightT> &att) {
  os << "(" << att.get_weight() << ")";
  return os;
}

template<typename WeightT>
std::istream& operator>>(std::istream &is, MinimalAttributes<WeightT> &att) {
  is >> att.weight_;
  if (!is) { att.weight_ = WeightTraits<WeightT>::one(); }
  return is;
}

template<typename WeightT>
MinimalAttributes<WeightT>::MinimalAttributes() 
  : weight_(WeightTraits<WeightT>::one()) {
}

template<typename WeightT>
MinimalAttributes<WeightT>::~MinimalAttributes() {}

template<typename WeightT>
MinimalAttributes<WeightT>::MinimalAttributes(const WeightT &w) 
  : weight_(w) {
}

template<typename WeightT>
inline const WeightT& MinimalAttributes<WeightT>::get_weight() const {
  return weight_;
}

template<typename WeightT>
inline void MinimalAttributes<WeightT>::set_weight(const WeightT &w) {
  weight_ = w;
}

// Attribute type functions
// ------------------------

// return a AttrT object that is the result of merging the arguments
template<typename AttrT> 
AttrT merge(const AttrT &lhs, const AttrT &rhs) {
  AttrT result = lhs; 
  result.merge(rhs);
  return result;
}

}

#endif
