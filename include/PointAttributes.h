#ifndef BBRCITKDE_POINTATTRIBUTES_H__
#define BBRCITKDE_POINTATTRIBUTES_H__

#include <iostream>
#include <KdeTraits.h>

namespace bbrcit {

template<typename WeightT> class MinimalAttributes;

template<typename WeightT> 
std::istream& operator>>(std::istream&, MinimalAttributes<WeightT>&);

template<typename WeightT> 
std::ostream& operator<<(std::ostream&, const MinimalAttributes<WeightT>&);

template<typename WeightT>
class MinimalAttributes {

  friend std::istream& operator>><>(std::istream&, MinimalAttributes<WeightT>&);

  public: 

    using WeightType = WeightT;

    MinimalAttributes();
    MinimalAttributes(const WeightT&);
    ~MinimalAttributes();

    const WeightT& get_weight() const;
    void set_weight(const WeightT&);

  private: 
    WeightT weight_;
};

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

}

#endif
