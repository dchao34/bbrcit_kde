#ifndef BBRCITKDE_DEFAULTNODEATTRIBUTES_H__
#define BBRCITKDE_DEFAULTNODEATTRIBUTES_H__

#include <iostream>

#include <AttributesInterface.h>

namespace bbrcit {

// DefaultNodeAttributes<>
// -----------------------

template<typename PointAttrT> class DefaultNodeAttributes;

// prints the attribute to ostream. 
template<typename PointAttrT> 
std::ostream& operator<<(std::ostream&, const DefaultNodeAttributes<PointAttrT>&);

template<typename PointAttrT> 
void swap(DefaultNodeAttributes<PointAttrT>&, DefaultNodeAttributes<PointAttrT>&);

// DefaultNodeAttributes<> contains a single attribute PointAttrT that is itself an attribute. 
// Example usage: In a Kdtree, PointAttrT is the attribute type of the data points. The natural
// default attribute for a node is simply the merger of all the leaves in its subtree. 
template<typename PointAttrT>
class DefaultNodeAttributes {

  public: 

    using PointAttributeType = PointAttrT;
    friend void swap<>(DefaultNodeAttributes<PointAttrT>&, DefaultNodeAttributes<PointAttrT>&);

    // required functions for the node attribute interface
    // + extract_point_attributes: takes a PointT argument representing a single 
    //   point in the Kdtree and returns an NodeAttributes<> object. 
    template<typename PointT> 
      static DefaultNodeAttributes extract_point_attributes(const PointT&);

    // default constructor: default initializes PointAttrT.
    DefaultNodeAttributes();

    // one argument constructor: initializes the attribute to a user defined value.  
    DefaultNodeAttributes(const PointAttrT&);

    // copy-control. 
    DefaultNodeAttributes(const DefaultNodeAttributes<PointAttrT>&) = default;
    DefaultNodeAttributes(DefaultNodeAttributes<PointAttrT>&&) = default;
    DefaultNodeAttributes<PointAttrT>& operator=(const DefaultNodeAttributes<PointAttrT>&) = default;
    DefaultNodeAttributes<PointAttrT>& operator=(DefaultNodeAttributes<PointAttrT>&&) = default;
    ~DefaultNodeAttributes();

    // access the attribute object
    const PointAttrT& get_point_attributes() const;
    void set_point_attributes(const PointAttrT&);

    // required functions for the attribute interface
    DefaultNodeAttributes<PointAttrT>& merge(const DefaultNodeAttributes<PointAttrT>&);

  private: 
    PointAttrT attributes_;
};

template<typename PointAttrT> 
  template<typename PointT> 
inline DefaultNodeAttributes<PointAttrT> 
DefaultNodeAttributes<PointAttrT>::extract_point_attributes(
    const PointT &p) {
  return DefaultNodeAttributes(p.get_attributes());
}

template<typename PointAttrT> 
inline void swap(DefaultNodeAttributes<PointAttrT> &lhs, DefaultNodeAttributes<PointAttrT> &rhs) {
  using std::swap; swap(lhs.attributes_, rhs.attributes_);
}

template<typename PointAttrT>
inline DefaultNodeAttributes<PointAttrT>& 
DefaultNodeAttributes<PointAttrT>::merge(const DefaultNodeAttributes<PointAttrT> &rhs) {
  attributes_.merge(rhs.attributes_);
  return *this;
}

template<typename PointAttrT> 
std::ostream& operator<<(std::ostream &os, const DefaultNodeAttributes<PointAttrT> &att) {
  os << att.get_point_attributes();
  return os;
}

template<typename PointAttrT>
DefaultNodeAttributes<PointAttrT>::DefaultNodeAttributes() 
  : attributes_() {
}

template<typename PointAttrT>
DefaultNodeAttributes<PointAttrT>::~DefaultNodeAttributes() {}

template<typename PointAttrT>
DefaultNodeAttributes<PointAttrT>::DefaultNodeAttributes(const PointAttrT &att) 
  : attributes_(att) {
}

template<typename PointAttrT>
inline const PointAttrT& DefaultNodeAttributes<PointAttrT>::get_point_attributes() const {
  return attributes_;
}

template<typename PointAttrT>
inline void DefaultNodeAttributes<PointAttrT>::set_point_attributes(const PointAttrT &att) {
  attributes_ = att;
}

}

#endif
