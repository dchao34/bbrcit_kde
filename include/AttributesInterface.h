#ifndef BBRCITKDE_ATTRIBUTESINTERFACE_H__
#define BBRCITKDE_ATTRIBUTESINTERFACE_H__

namespace bbrcit {

// Interface: All Attributes<> types (AttrT) must meet the following requirements:
// + Member functions:
//     + AttrT& merge(const AttrT &rhs): Defines how to merge rhs into *this. 
//       returns *this. 

// return a AttrT object that is the result of merging the arguments
template<typename AttrT> 
AttrT merge(const AttrT &lhs, const AttrT &rhs) {
  AttrT result = lhs; 
  result.merge(rhs);
  return result;
}

}

#endif
