#ifndef BBRCITKDE_EMPTYNODEATTRIBUTES_H__
#define BBRCITKDE_EMPTYNODEATTRIBUTES_H__

#include <iostream>

#include <AttributesInterface.h>

namespace bbrcit {

// EmptyNodeAttributes<>
// ---------------------

template<typename PointAttrT> class EmptyNodeAttributes;

// prints the attribute to ostream. 
template<typename PointAttrT> 
inline std::ostream& operator<<(std::ostream &os, const EmptyNodeAttributes<PointAttrT>&) { os << "()"; return os; }

// EmptyNodeAttributes<> represents the empty attribute. 
template<typename PointAttrT>
class EmptyNodeAttributes {

  public: 

    using PointAttributeType = PointAttrT;

    // default constructor
    EmptyNodeAttributes() = default;
    EmptyNodeAttributes(const EmptyNodeAttributes<PointAttrT>&) = default;
    EmptyNodeAttributes(EmptyNodeAttributes<PointAttrT>&&) = default;
    EmptyNodeAttributes<PointAttrT>& operator=(const EmptyNodeAttributes<PointAttrT>&) = default;
    EmptyNodeAttributes<PointAttrT>& operator=(EmptyNodeAttributes<PointAttrT>&&) = default;
    ~EmptyNodeAttributes() = default;

    // required functions for the attribute interface
    EmptyNodeAttributes<PointAttrT>& merge(const EmptyNodeAttributes<PointAttrT>&) { return *this; }

    // required functions for the node attribute interface
    template<typename PointT> void configure_for_leaf(const PointT&) { return; }

  private: 
    PointAttrT attributes_;
};

}

#endif
