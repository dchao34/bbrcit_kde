#ifndef BBRCITKDE_KDETRAITS_H__
#define BBRCITKDE_KDETRAITS_H__

#ifdef __CUDACC__
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_CALLABLE
#endif

#include <Point2d.h>

namespace bbrcit {

template<typename T> class ConstantTraits;

template<> 
class ConstantTraits<double> {
  public: 
    using Type = double;
    CUDA_CALLABLE static Type zero() { return 0.0; } 
    CUDA_CALLABLE static Type one() { return 1.0; } 
};

template<> 
class ConstantTraits<float> {
  public: 
    using Type = float;
    CUDA_CALLABLE static Type zero() { return 0.0f; } 
    CUDA_CALLABLE static Type one() { return 1.0f; } 
};

template<> 
class ConstantTraits<int> {
  public: 
    using Type = int;
    CUDA_CALLABLE static Type zero() { return 0; } 
    CUDA_CALLABLE static Type one() { return 1; } 
};

template<> 
class ConstantTraits<unsigned> {
  public: 
    using Type = unsigned;
    CUDA_CALLABLE static Type zero() { return 0; } 
    CUDA_CALLABLE static Type one() { return 1; } 
};

#ifdef __CUDACC__

template<int D, typename FloatT> class DevicePointTraits;

template<> 
class DevicePointTraits<2, float> {
  public: 
    using Type = Point2d<float>;
};

template<> 
class DevicePointTraits<2, double> {
  public: 
    using Type = Point2d<double>;
};

#endif


}

#endif
