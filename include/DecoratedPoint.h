#ifndef BBRCITKDE_DECORATEDPOINT_H__
#define BBRCITKDE_DECORATEDPOINT_H__

#include <vector>
#include <utility>
#include <exception>
#include <initializer_list>
#include <iostream>

#include <Point.h>
#include <PointAttributes.h>

namespace bbrcit {

// API
// ---

template <int Dim, typename AttrT, typename FloatT> class DecoratedPoint;

// prints a DecoratedPoint<> as { (1), (2) } to os, where (1) is operator<< of
// the PointT and (2) is operator<< of AttributesT.
template <int Dim, typename AttrT, typename FloatT>
std::ostream& operator<<(std::ostream&, const DecoratedPoint<Dim,AttrT,FloatT>&);

// reads in two steps: (1) operator>> of PointT. (2) operator>> of AttributesT. 
template <int Dim, typename AttrT, typename FloatT>
std::istream& operator>>(std::istream&, DecoratedPoint<Dim,AttrT,FloatT>&);

template <int Dim, typename AttrT, typename FloatT>
void swap(DecoratedPoint<Dim,AttrT,FloatT>&, DecoratedPoint<Dim,AttrT,FloatT>&);

// DecoratedPoint<> models a weighted point in Dim-dimensional space. 
template <int Dim, typename AttrT=MinimalAttributes<int>, typename FloatT=double> 
class DecoratedPoint {

  public:

    using PointT = Point<Dim, FloatT>;
    using AttributesT = AttrT;

    friend std::istream& operator>><>(std::istream&, DecoratedPoint<Dim,AttrT,FloatT>&);
    friend std::ostream& operator<<<>(std::ostream&, const DecoratedPoint<Dim,AttrT,FloatT>&);
    friend void swap<>(DecoratedPoint<Dim,AttrT,FloatT>&, DecoratedPoint<Dim,AttrT,FloatT>&);

    // default constructor: creates point at the origin with default attributes.
    DecoratedPoint();

    // PointT constructor: creates DecoratedPoint<> centered at PointT. 
    // attribute is configured as follows:
    // (1-2) Default attributes: AttributesT().
    // (3) User specifies the attribute in the argument. 
    DecoratedPoint(const PointT&);
    DecoratedPoint(PointT&&);
    DecoratedPoint(const PointT&, const AttrT&);

    // copy-control. operator= uses copy and swap. 
    DecoratedPoint(const DecoratedPoint<Dim,AttrT,FloatT>&);
    DecoratedPoint(DecoratedPoint<Dim,AttrT,FloatT>&&);
    DecoratedPoint<Dim,AttrT,FloatT>& operator=(DecoratedPoint<Dim,AttrT,FloatT>);
    ~DecoratedPoint();

    // returns the attribute or point
    const AttrT& get_attributes() const;
    const PointT& get_point() const;

    // set the weight or point
    void set_attributes(const AttrT&);
    void set_point(const PointT&);

    // returns the ith coordinate of the point. 
    const FloatT& operator[](size_t i) const;
    FloatT& operator[](size_t);

  private: 

    // represents this DecoratedPoint<>'s location in Dim-dimensional euclidean space
    PointT point_;

    // attributes of this DecoratedPoint<> object. 
    AttrT attr_;
};

// Implementations
// ---------------

template <int Dim, typename AttrT, typename FloatT>
void swap(DecoratedPoint<Dim,AttrT,FloatT> &lhs, DecoratedPoint<Dim,AttrT,FloatT> &rhs) {
  using std::swap;
  swap(lhs.point_, rhs.point_);
  swap(lhs.attr_, rhs.attr_);
}

template <int Dim, typename AttrT, typename FloatT>
inline const FloatT& 
DecoratedPoint<Dim,AttrT,FloatT>::operator[](size_t i) const {
  return point_[i];
}

template <int Dim, typename AttrT, typename FloatT>
inline FloatT& DecoratedPoint<Dim,AttrT,FloatT>::operator[](size_t i) {
  return const_cast<FloatT&>(
      static_cast<const DecoratedPoint<Dim,AttrT,FloatT>&>(*this)[i]
  );
}

template <int Dim, typename AttrT, typename FloatT>
std::ostream& operator<<(
    std::ostream &os, 
    const DecoratedPoint<Dim,AttrT,FloatT> &p) {
  os << "{ " << p.point_ << ", " << p.attr_ << " }";
  return os;
}

template <int Dim, typename AttrT, typename FloatT>
std::istream& operator>>(
    std::istream &is, 
    DecoratedPoint<Dim,AttrT,FloatT> &p) {
  is >> p.point_; is >> p.attr_;
  if (!is) { p = DecoratedPoint<Dim,AttrT,FloatT>(); }
  return is;
}

template <int Dim, typename AttrT, typename FloatT> 
inline const AttrT& 
DecoratedPoint<Dim,AttrT,FloatT>::get_attributes() const 
{ return attr_; }

template <int Dim, typename AttrT, typename FloatT> 
inline const typename DecoratedPoint<Dim,AttrT,FloatT>::PointT& 
DecoratedPoint<Dim,AttrT,FloatT>::get_point() const 
{ return point_; }

template <int Dim, typename AttrT, typename FloatT> 
inline void 
DecoratedPoint<Dim,AttrT,FloatT>::set_attributes(const AttrT &a)
{ attr_ = a; }

template <int Dim, typename AttrT, typename FloatT> 
inline void 
DecoratedPoint<Dim,AttrT,FloatT>::set_point(
  const typename DecoratedPoint<Dim,AttrT,FloatT>::PointT &p) 
{ point_ = p; }

template <int Dim, typename AttrT, typename FloatT> 
DecoratedPoint<Dim, AttrT, FloatT>::DecoratedPoint() 
  : point_(), attr_() {
}

template <int Dim, typename AttrT, typename FloatT> 
DecoratedPoint<Dim, AttrT, FloatT>::~DecoratedPoint() {}

template <int Dim, typename AttrT, typename FloatT> 
DecoratedPoint<Dim, AttrT, FloatT>::DecoratedPoint(const PointT &p) 
  : point_(p), attr_() {
}

template <int Dim, typename AttrT, typename FloatT> 
DecoratedPoint<Dim, AttrT, FloatT>::DecoratedPoint(PointT &&p) 
  : point_(std::move(p)), attr_() {
}

template <int Dim, typename AttrT, typename FloatT> 
DecoratedPoint<Dim, AttrT, FloatT>::DecoratedPoint(const PointT &p, const AttrT &a) 
  : point_(p), attr_(a) {
}

template <int Dim, typename AttrT, typename FloatT> 
DecoratedPoint<Dim, AttrT, FloatT>::DecoratedPoint(const DecoratedPoint<Dim,AttrT,FloatT> &p) 
  : point_(p.point_), attr_(p.attr_) {
}

template <int Dim, typename AttrT, typename FloatT> 
DecoratedPoint<Dim, AttrT, FloatT>::DecoratedPoint(DecoratedPoint<Dim,AttrT,FloatT> &&p) 
  : point_(std::move(p.point_)), attr_(std::move(p.attr_)) {
}

template <int Dim, typename AttrT, typename FloatT> 
inline DecoratedPoint<Dim,AttrT,FloatT>& 
DecoratedPoint<Dim,AttrT,FloatT>::operator=(DecoratedPoint<Dim,AttrT,FloatT> rhs) {
  swap(*this, rhs); return *this;
}


}
#endif
