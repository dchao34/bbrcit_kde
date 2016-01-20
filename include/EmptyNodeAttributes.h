#ifndef BBRCITKDE_EMPTYNODEATTRIBUTES_H__
#define BBRCITKDE_EMPTYNODEATTRIBUTES_H__

#include <iostream>

#include <AttributesInterface.h>

namespace bbrcit {

// EmptyNodeAttributes
// ---------------------

class EmptyNodeAttributes;

// prints the attribute to ostream. 
inline std::ostream& operator<<(std::ostream &os, const EmptyNodeAttributes&) { os << "()"; return os; }

// EmptyNodeAttributes represents the empty attribute. 
class EmptyNodeAttributes {

  public: 

    // required function for the node attribute interface
    template<typename PointT> 
      static EmptyNodeAttributes extract_point_attributes(const PointT&);

    // default constructor
    EmptyNodeAttributes() = default;
    EmptyNodeAttributes(const EmptyNodeAttributes&) = default;
    EmptyNodeAttributes(EmptyNodeAttributes&&) = default;
    EmptyNodeAttributes& operator=(const EmptyNodeAttributes&) = default;
    EmptyNodeAttributes& operator=(EmptyNodeAttributes&&) = default;
    ~EmptyNodeAttributes() = default;

    // required functions for the attribute interface
    EmptyNodeAttributes& merge(const EmptyNodeAttributes&) { return *this; }

};

template<typename PointT>
inline EmptyNodeAttributes 
EmptyNodeAttributes::extract_point_attributes(const PointT&) 
{ return EmptyNodeAttributes(); }

}

#endif
