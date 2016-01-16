#ifndef BBRCITKDE_DECORATEDPOINT_H__
#define BBRCITKDE_DECORATEDPOINT_H__

#include <vector>
#include <utility>
#include <exception>
#include <initializer_list>
#include <iostream>

#include <Point.h>

namespace bbrcit {

// API
// ---

template <int Dim, typename AttrT> class DecoratedPoint;

// prints a DecoratedPoint<> as { (1), (2) } to os, where (1) is operator<< of
// the PointT and (2) is operator<< of AttributesT.
template <int Dim, typename AttrT>
std::ostream& operator<<(std::ostream&, const DecoratedPoint<Dim,AttrT>&);

// reads in two steps: (1) operator>> of PointT. (2) operator>> of AttributesT. 
template <int Dim, typename AttrT>
std::istream& operator>>(std::istream&, DecoratedPoint<Dim,AttrT>&);

template <int Dim, typename AttrT>
void swap(DecoratedPoint<Dim,AttrT>&, DecoratedPoint<Dim,AttrT>&);

// DecoratedPoint<> models a point in Dim-dimensional space with additional attributes attached.
template <int Dim, typename AttrT> 
class DecoratedPoint {

  public:

    using FloatType = typename AttrT::FloatType;
    using PointT = Point<Dim, FloatType>;
    using AttributesT = AttrT;

    friend std::istream& operator>><>(std::istream&, DecoratedPoint<Dim,AttrT>&);
    friend std::ostream& operator<<<>(std::ostream&, const DecoratedPoint<Dim,AttrT>&);
    friend void swap<>(DecoratedPoint<Dim,AttrT>&, DecoratedPoint<Dim,AttrT>&);

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
    DecoratedPoint(const DecoratedPoint<Dim,AttrT>&);
    DecoratedPoint(DecoratedPoint<Dim,AttrT>&&);
    DecoratedPoint<Dim,AttrT>& operator=(DecoratedPoint<Dim,AttrT>);
    ~DecoratedPoint();

    // returns the attribute or point
    const AttrT& get_attributes() const;
    const PointT& get_point() const;

    // set the attributes or point
    void set_attributes(const AttrT&);
    void set_point(const PointT&);

    // returns the ith coordinate of the point. 
    const FloatType& operator[](size_t i) const;
    FloatType& operator[](size_t);

  private: 

    // represents this DecoratedPoint<>'s location in Dim-dimensional euclidean space
    PointT point_;

    // attributes of this DecoratedPoint<> object. 
    AttrT attr_;
};

// Implementations
// ---------------

template <int Dim, typename AttrT>
void swap(DecoratedPoint<Dim,AttrT> &lhs, DecoratedPoint<Dim,AttrT> &rhs) {
  using std::swap;
  swap(lhs.point_, rhs.point_);
  swap(lhs.attr_, rhs.attr_);
}

template <int Dim, typename AttrT>
inline const typename DecoratedPoint<Dim,AttrT>::FloatType& 
DecoratedPoint<Dim,AttrT>::operator[](size_t i) const {
  return point_[i];
}

template <int Dim, typename AttrT>
inline typename DecoratedPoint<Dim,AttrT>::FloatType& 
DecoratedPoint<Dim,AttrT>::operator[](size_t i) {
  return const_cast<FloatType&>(
      static_cast<const DecoratedPoint<Dim,AttrT>&>(*this)[i]
  );
}

template <int Dim, typename AttrT>
std::ostream& operator<<(
    std::ostream &os, 
    const DecoratedPoint<Dim,AttrT> &p) {
  os << "{ " << p.point_ << ", " << p.attr_ << " }";
  return os;
}

template <int Dim, typename AttrT>
std::istream& operator>>(
    std::istream &is, 
    DecoratedPoint<Dim,AttrT> &p) {
  is >> p.point_; is >> p.attr_;
  if (!is) { p = DecoratedPoint<Dim,AttrT>(); }
  return is;
}

template <int Dim, typename AttrT> 
inline const AttrT& 
DecoratedPoint<Dim,AttrT>::get_attributes() const 
{ return attr_; }

template <int Dim, typename AttrT> 
inline const typename DecoratedPoint<Dim,AttrT>::PointT& 
DecoratedPoint<Dim,AttrT>::get_point() const 
{ return point_; }

template <int Dim, typename AttrT> 
inline void 
DecoratedPoint<Dim,AttrT>::set_attributes(const AttrT &a)
{ attr_ = a; }

template <int Dim, typename AttrT> 
inline void 
DecoratedPoint<Dim,AttrT>::set_point(
  const typename DecoratedPoint<Dim,AttrT>::PointT &p) 
{ point_ = p; }

template <int Dim, typename AttrT> 
DecoratedPoint<Dim, AttrT>::DecoratedPoint() 
  : point_(), attr_() {
}

template <int Dim, typename AttrT> 
DecoratedPoint<Dim, AttrT>::~DecoratedPoint() {}

template <int Dim, typename AttrT> 
DecoratedPoint<Dim, AttrT>::DecoratedPoint(const PointT &p) 
  : point_(p), attr_() {
}

template <int Dim, typename AttrT> 
DecoratedPoint<Dim, AttrT>::DecoratedPoint(PointT &&p) 
  : point_(std::move(p)), attr_() {
}

template <int Dim, typename AttrT> 
DecoratedPoint<Dim, AttrT>::DecoratedPoint(const PointT &p, const AttrT &a) 
  : point_(p), attr_(a) {
}

template <int Dim, typename AttrT> 
DecoratedPoint<Dim, AttrT>::DecoratedPoint(const DecoratedPoint<Dim,AttrT> &p) 
  : point_(p.point_), attr_(p.attr_) {
}

template <int Dim, typename AttrT> 
DecoratedPoint<Dim, AttrT>::DecoratedPoint(DecoratedPoint<Dim,AttrT> &&p) 
  : point_(std::move(p.point_)), attr_(std::move(p.attr_)) {
}

template <int Dim, typename AttrT> 
inline DecoratedPoint<Dim,AttrT>& 
DecoratedPoint<Dim,AttrT>::operator=(DecoratedPoint<Dim,AttrT> rhs) {
  swap(*this, rhs); return *this;
}


}
#endif
