#ifndef BBRCITKDE_POINTW_H__
#define BBRCITKDE_POINTW_H__

#include <vector>
#include <utility>
#include <exception>
#include <initializer_list>
#include <iostream>

#include <Point.h>
#include <KdeTraits.h>

namespace bbrcit {

// API
// ---

template <int Dim, typename WeightT, typename FloatT> class PointW;

// prints a PointW<> as { (coord1,...,coordDim), weight } to os. 
template <int Dim, typename WeightT, typename FloatT>
std::ostream& operator<<(std::ostream&, const PointW<Dim,WeightT,FloatT>&);

// reads D+1 whitespace separated values such that the first D initialize
// the point coordiates and the last initializes the weight.
template <int Dim, typename WeightT, typename FloatT>
std::istream& operator>>(std::istream&, PointW<Dim,WeightT,FloatT>&);

template <int Dim, typename WeightT, typename FloatT>
void swap(PointW<Dim,WeightT,FloatT>&, PointW<Dim,WeightT,FloatT>&);

// PointW<> models a weighted point in Dim-dimensional space. 
template <int Dim, typename WeightT, typename FloatT=double> 
class PointW {

  public:

    using PointT = Point<Dim, FloatT>;

    friend std::istream& operator>><>(std::istream&, PointW<Dim,WeightT,FloatT>&);
    friend void swap<>(PointW<Dim,WeightT,FloatT>&, PointW<Dim,WeightT,FloatT>&);

    // default constructor: creates point at the origin with weight 1.
    PointW();

    // PointT constructor: creates PointW<> centered at PointT. The weight 
    // is configured as follows:
    // (1-2) WeightTraits<WeightT>::one().
    // (3) User specifies the weight in the argument. 
    PointW(const PointT&);
    PointW(PointT&&);
    PointW(const PointT&, const WeightT&);

    // copy-control. operator= uses copy and swap. 
    PointW(const PointW<Dim,WeightT,FloatT>&);
    PointW(PointW<Dim,WeightT,FloatT>&&);
    PointW<Dim,WeightT,FloatT>& operator=(PointW<Dim,WeightT,FloatT>);
    ~PointW();

    // returns the weight or point
    const WeightT& get_weight() const;
    const PointT& get_point() const;

    // set the weight or point
    void set_weight(const WeightT&);
    void set_point(const PointT&);

    // returns the ith coordinate of the point. 
    const FloatT& operator[](size_t i) const;
    FloatT& operator[](size_t);

  private: 

    // represents this PointW<>'s location in Dim-dimensional euclidean space
    PointT point_;

    // weight of this PointW<> object. 
    WeightT weight_;
};

// Implementations
// ---------------

template <int Dim, typename WeightT, typename FloatT>
void swap(PointW<Dim,WeightT,FloatT> &lhs, PointW<Dim,WeightT,FloatT> &rhs) {
  using std::swap;
  swap(lhs.point_, rhs.point_);
  swap(lhs.weight_, rhs.weight_);
}

template <int Dim, typename WeightT, typename FloatT>
inline const FloatT& 
PointW<Dim,WeightT,FloatT>::operator[](size_t i) const {
  return point_[i];
}

template <int Dim, typename WeightT, typename FloatT>
inline FloatT& PointW<Dim,WeightT,FloatT>::operator[](size_t i) {
  return const_cast<FloatT&>(
      static_cast<const PointW<Dim,WeightT,FloatT>&>(*this)[i]
  );
}

template <int Dim, typename WeightT, typename FloatT>
std::ostream& operator<<(
    std::ostream &os, 
    const PointW<Dim,WeightT,FloatT> &p) {
  os << "{ " << p.get_point() << ", " << p.get_weight() << " }";
  return os;
}

template <int Dim, typename WeightT, typename FloatT>
std::istream& operator>>(
    std::istream &is, 
    PointW<Dim,WeightT,FloatT> &p) {
  is >> p.point_; is >> p.weight_;
  if (!is) { p = PointW<Dim,WeightT,FloatT>(); }
  return is;
}

template <int Dim, typename WeightT, typename FloatT>
std::istream& operator>>(std::istream&, PointW<Dim,WeightT,FloatT>&);

template <int Dim, typename WeightT, typename FloatT> 
inline const WeightT& 
PointW<Dim,WeightT,FloatT>::get_weight() const 
{ return weight_; }

template <int Dim, typename WeightT, typename FloatT> 
inline const typename PointW<Dim,WeightT,FloatT>::PointT& 
PointW<Dim,WeightT,FloatT>::get_point() const 
{ return point_; }

template <int Dim, typename WeightT, typename FloatT> 
inline void 
PointW<Dim,WeightT,FloatT>::set_weight(const WeightT &w)
{ weight_ = w; }

template <int Dim, typename WeightT, typename FloatT> 
inline void 
PointW<Dim,WeightT,FloatT>::set_point(
  const typename PointW<Dim,WeightT,FloatT>::PointT &p) 
{ point_ = p; }

template <int Dim, typename WeightT, typename FloatT> 
PointW<Dim, WeightT, FloatT>::PointW() 
  : point_(), weight_(WeightTraits<WeightT>::one()) {
}

template <int Dim, typename WeightT, typename FloatT> 
PointW<Dim, WeightT, FloatT>::~PointW() {}

template <int Dim, typename WeightT, typename FloatT> 
PointW<Dim, WeightT, FloatT>::PointW(const PointT &p) 
  : point_(p), weight_(WeightTraits<WeightT>::one()) {
}

template <int Dim, typename WeightT, typename FloatT> 
PointW<Dim, WeightT, FloatT>::PointW(PointT &&p) 
  : point_(std::move(p)), weight_(WeightTraits<WeightT>::one()) {
}

template <int Dim, typename WeightT, typename FloatT> 
PointW<Dim, WeightT, FloatT>::PointW(const PointT &p, const WeightT &w) 
  : point_(p), weight_(w) {
}

template <int Dim, typename WeightT, typename FloatT> 
PointW<Dim, WeightT, FloatT>::PointW(const PointW<Dim,WeightT,FloatT> &p) 
  : point_(p.point_), weight_(p.weight_) {
}

template <int Dim, typename WeightT, typename FloatT> 
PointW<Dim, WeightT, FloatT>::PointW(PointW<Dim,WeightT,FloatT> &&p) 
  : point_(std::move(p.point_)), weight_(std::move(p.weight_)) {
}

template <int Dim, typename WeightT, typename FloatT> 
inline PointW<Dim,WeightT,FloatT>& 
PointW<Dim,WeightT,FloatT>::operator=(PointW<Dim,WeightT,FloatT> rhs) {
  swap(*this, rhs); return *this;
}


}
#endif
