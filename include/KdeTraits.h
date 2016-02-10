#ifndef BBRCITKDE_KDETRAITS_H__
#define BBRCITKDE_KDETRAITS_H__

namespace bbrcit {

template<typename T> class ConstantTraits;

template<> 
class ConstantTraits<double> {
  public: 
    using Type = double;
    static Type zero() { return 0.0; } 
    static Type one() { return 1.0; } 
};

template<> 
class ConstantTraits<float> {
  public: 
    using Type = float;
    static Type zero() { return 0.0; } 
    static Type one() { return 1.0; } 
};

template<> 
class ConstantTraits<int> {
  public: 
    using Type = int;
    static Type zero() { return 0; } 
    static Type one() { return 1; } 
};

template<> 
class ConstantTraits<unsigned> {
  public: 
    using Type = unsigned;
    static Type zero() { return 0; } 
    static Type one() { return 1; } 
};

}

#endif
