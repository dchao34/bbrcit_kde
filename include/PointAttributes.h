#ifndef BBRCITKDE_POINTATTRIBUTES_H__
#define BBRCITKDE_POINTATTRIBUTES_H__

#include <iostream>
#include <KdeTraits.h>

namespace bbrcit {

template<typename WeightT, typename FloatT> class MinimalAttributes;

// prints MinimalAttributes<> as { (weight) } to os. 
template<typename WeightT, typename FloatT> 
std::istream& operator>>(std::istream&, MinimalAttributes<WeightT,FloatT>&);

// reads in a single token
template<typename WeightT, typename FloatT> 
std::ostream& operator<<(std::ostream&, const MinimalAttributes<WeightT,FloatT>&);

// MinimalAttributes<> contains the minimal attributes required of a Kdtree data point. 
template<typename WeightT, typename FloatT>
class MinimalAttributes {

  friend std::istream& operator>><>(std::istream&, MinimalAttributes<WeightT,FloatT>&);

  public: 

    using WeightType = WeightT;
    using FloatType = FloatT;

    // default constructor: weight equal to 1.
    MinimalAttributes();

    // one-argument constructor: user configured weight. 
    MinimalAttributes(const WeightT&);
    ~MinimalAttributes();

    const WeightT& get_weight() const;
    void set_weight(const WeightT&);

  private: 
    WeightT weight_;
};

template<typename WeightT, typename FloatT> 
std::ostream& operator<<(std::ostream &os, const MinimalAttributes<WeightT,FloatT> &att) {
  os << "(" << att.get_weight() << ")";
  return os;
}

template<typename WeightT, typename FloatT>
std::istream& operator>>(std::istream &is, MinimalAttributes<WeightT,FloatT> &att) {
  is >> att.weight_;
  if (!is) { att.weight_ = WeightTraits<WeightT>::one(); }
  return is;
}

template<typename WeightT, typename FloatT>
MinimalAttributes<WeightT,FloatT>::MinimalAttributes() 
  : weight_(WeightTraits<WeightT>::one()) {
}

template<typename WeightT, typename FloatT>
MinimalAttributes<WeightT,FloatT>::~MinimalAttributes() {}

template<typename WeightT, typename FloatT>
MinimalAttributes<WeightT,FloatT>::MinimalAttributes(const WeightT &w) 
  : weight_(w) {
}

template<typename WeightT, typename FloatT>
inline const WeightT& MinimalAttributes<WeightT,FloatT>::get_weight() const {
  return weight_;
}

template<typename WeightT, typename FloatT>
inline void MinimalAttributes<WeightT,FloatT>::set_weight(const WeightT &w) {
  weight_ = w;
}

}

#endif
