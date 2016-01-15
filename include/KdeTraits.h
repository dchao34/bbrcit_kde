#ifndef BBRCITKDE_KDETRAITS_H__
#define BBRCITKDE_KDETRAITS_H__

namespace bbrcit {

template<typename T> class WeightTraits;

template<> 
class WeightTraits<double> {
  public: 
    using WeightT = double;
    static WeightT zero() { return 0; } 
    static WeightT one() { return 1; } 
};

template<> 
class WeightTraits<float> {
  public: 
    using WeightT = float;
    static WeightT zero() { return 0; } 
    static WeightT one() { return 1; } 
};

template<> 
class WeightTraits<int> {
  public: 
    using WeightT = int;
    static WeightT zero() { return 0; } 
    static WeightT one() { return 1; } 
};

template<> 
class WeightTraits<unsigned> {
  public: 
    using WeightT = unsigned;
    static WeightT zero() { return 0; } 
    static WeightT one() { return 1; } 
};

}

#endif
