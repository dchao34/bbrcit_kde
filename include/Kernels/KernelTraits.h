#ifndef BBRCITKDE_KERNELTRAITS_H__
#define BBRCITKDE_KERNELTRAITS_H__

#ifdef __CUDACC__
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_CALLABLE
#endif

#define _USE_MATH_DEFINES

#include <cmath>

namespace bbrcit {

// Epanechnikov
// ------------

template<int D, typename T> 
class EpanechnikovTraits {
  public:
    static T normalization(T h) { 
      return 0.5 * (D+2) / (std::pow(M_PI, D/2.0) / std::tgamma(1+D/2.0))
                     / std::pow(h, D);
    }

    // evaluates max(1-point_arg, 0)
    static T kernel(T point_arg) {
      return std::max(1.0 - point_arg, 0.0);
    }
};

template<> 
class EpanechnikovTraits<1,float> {

  public:

    // 0.5 * 3 / (sqrt(pi) * gamma(3/2)) / h = 0.75 / h
    CUDA_CALLABLE static float normalization(float h) { 
      return 0.75f / h;
    }

    CUDA_CALLABLE static float kernel(float point_arg) {
      return fmaxf(1.0f - point_arg, 0.0f);
    }
};

template<> 
class EpanechnikovTraits<1,double> {
  
  public:

    CUDA_CALLABLE static double normalization(double h) { 
      return 0.75 / h;
    }

    CUDA_CALLABLE static double kernel(double point_arg) {
      return fmax(1.0 - point_arg, 0.0);
    }
};

template<> 
class EpanechnikovTraits<2,float> {

  public:

    // 0.5 * 4 / (pi*h*h) = 2 / (pi*h*h) =  0.63661977236758134307553 / (h*h)
    CUDA_CALLABLE static float normalization(float h) { 
      return 0.6366197723f / (h*h);
    }

    CUDA_CALLABLE static float kernel(float point_arg) {
      return fmaxf(1.0f - point_arg, 0.0f);
    }
};

template<> 
class EpanechnikovTraits<2,double> {
  public:

    CUDA_CALLABLE static double normalization(double h) { 
      return 0.6366197723675813431 / (h*h);
    }

    CUDA_CALLABLE static double kernel(double point_arg) {
      return fmax(1.0 - point_arg, 0.0);
    }
};





// Gaussian
// --------

template<int D, typename T> 
class GaussianTraits {
  public:
    static T normalization(T h) { 
      return pow(2*M_PI, -D/2.0) / pow(h, D);
    }

    // evaluates exp(-0.5 * point_arg)
    static T kernel(T point_arg) {
      return exp(-0.5 * point_arg);
    }
};

template<>
class GaussianTraits<1,float> {
  public:

    // evaluates 1 / (sqrt(2pi)) / h; 1 / sqrt(2pi) = 0.39894228040143267794
    CUDA_CALLABLE static float normalization(float h) { 
      return 0.3989422804f / h;
    }

    CUDA_CALLABLE static float kernel(float point_arg) {
      return expf(-0.5f * point_arg);
    }
};

template<>
class GaussianTraits<2,float> {
  public:

    // evaluates 1 / (2pi*h*h); 1 / (2pi) = 0.15915494309189533577
    CUDA_CALLABLE static float normalization(float h) { 
      return 0.1591549431f / (h*h);
    }

    CUDA_CALLABLE static float kernel(float point_arg) {
      return expf(-0.5f * point_arg);
    }
};

template<>
class GaussianTraits<1,double> {
  public:

    // evaluates 1 / (sqrt(2pi)) / h; 1 / sqrt(2pi) = 0.39894228040143267794
    CUDA_CALLABLE static double normalization(float h) { 
      return 0.3989422804014326779 / h;
    }

    CUDA_CALLABLE static double kernel(double point_arg) {
      return exp(-0.5 * point_arg);
    }
};

template<>
class GaussianTraits<2,double> {
  public:

    // evaluates 1 / (2pi*h*h); 1 / (2pi) = 0.15915494309189533577
    CUDA_CALLABLE static double normalization(double h) { 
      return 0.1591549430918953357 / (h*h);
    }

    CUDA_CALLABLE static double kernel(double point_arg) {
      return exp(-0.5 * point_arg);
    }
};


// GaussianProduct2d
// -----------------

template<typename T> 
class GaussianProduct2dTraits {

  public:
    // evaluates 1 / (2pi*hx*hy); 1 / (2pi) = 0.15915494309189533577
    static T normalization(T hx, T hy) { 
      return 0.1591549430918953357 / (hx*hy);
    }

    // evaluates exp(-0.5 * point_arg)
    static T kernel(T point_arg) {
      return exp(-0.5 * point_arg);
    }
};

template<> 
class GaussianProduct2dTraits<float> {

  public:
    // evaluates 1 / (2pi*hx*hy); 1 / (2pi) = 0.15915494309189533577
    CUDA_CALLABLE static float normalization(float hx, float hy) { 
      return 0.1591549431f / (hx*hy);
    }

    // evaluates exp(-0.5 * point_arg)
    CUDA_CALLABLE static float kernel(float point_arg) {
      return expf(-0.5f * point_arg);
    }
};

template<> 
class GaussianProduct2dTraits<double> {

  public:
    // evaluates 1 / (2pi*hx*hy); 1 / (2pi) = 0.15915494309189533577
    CUDA_CALLABLE static double normalization(double hx, double hy) { 
      return 0.1591549430918953357 / (hx*hy);
    }

    // evaluates exp(-0.5 * point_arg)
    CUDA_CALLABLE static double kernel(double point_arg) {
      return exp(-0.5 * point_arg);
    }
};

// EpanechnikovProduct2d
// ---------------------

template<typename T> 
class EpanechnikovProduct2dTraits {

  public:
    // evaluates (9/16) / (hx*hy); 9/16 = 0.5625
    static T normalization(T hx, T hy) { 
      return 0.5625 / (hx*hy);
    }

    // evaluates max(1-point_arg, 0)
    static T kernel(T point_arg) {
      return std::max(1.0 - point_arg, 0.0);
    }
};

template<> 
class EpanechnikovProduct2dTraits<float> {

  public:
    CUDA_CALLABLE static float normalization(float hx, float hy) { 
      return 0.5625f / (hx*hy);
    }

    CUDA_CALLABLE static float kernel(float point_arg) {
      return fmaxf(1.0f - point_arg, 0.0f);
    }
};

template<> 
class EpanechnikovProduct2dTraits<double> {

  public:
    CUDA_CALLABLE static double normalization(double hx, double hy) { 
      return 0.5625 / (hx*hy);
    }

    CUDA_CALLABLE static double kernel(double point_arg) {
      return fmax(1.0 - point_arg, 0.0);
    }
};


// EpanechnikovConvKernel1d
// ------------------------

template<typename T> 
class EpanechnikovConv1dTraits {

  public:
    // evaluates (3/160) / hx; 3/160 = 0.01875
    static T normalization(T hx) { 
      return 0.6 / hx;
    }

    // evaluates (2-x)^3*(x^2+6x+4) if 0<=|x|<2 and 0 otherwise. 
    static T kernel(T arg) {
      arg = std::abs(arg);
      T cubic = std::max(2-arg, 0);
      return 0.03125*cubic*cubic*cubic*(arg*arg+6*arg+4);
    }
};

template<> 
class EpanechnikovConv1dTraits<float> {

  public:
    CUDA_CALLABLE static float normalization(float hx) { 
      return 0.6f / hx;
    }

    CUDA_CALLABLE static float kernel(float arg) {
      arg = fabsf(arg);
      float cubic = fmaxf(2.0f-arg, 0.0f);
      return 0.03125f*cubic*cubic*cubic*(arg*arg+6.0f*arg+4.0f);
    }
};

template<> 
class EpanechnikovConv1dTraits<double> {

  public:
    CUDA_CALLABLE static double normalization(double hx) { 
      return 0.6 / hx;
    }

    CUDA_CALLABLE static double kernel(double arg) {
      arg = fabs(arg);
      double cubic = fmax(2.0-arg, 0.0);
      return 0.03125*cubic*cubic*cubic*(arg*arg+6.0*arg+4.0);
    }
};

// EpanechnikovProductConvKernel2d
// -------------------------------

template<typename T> 
class EpanechnikovProductConv2dTraits {

  public:
    // evaluates (3*3/(160*160)) / (hx*hy); 3*3/(160*160) = 0.0003515625
    static T normalization(T hx, T hy) { 
      return 0.36 / (hx*hy);
    }

    // evaluates (2-x)^3*(x^2+6x+4) if 0<=|x|<2 and 0 otherwise. 
    static T kernel(T arg) {
      arg = std::abs(arg);
      T cubic = std::max(2-arg, 0);
      return 0.03125*cubic*cubic*cubic*(arg*arg+6*arg+4);
    }
};

template<> 
class EpanechnikovProductConv2dTraits<float> {

  public:
    CUDA_CALLABLE static float normalization(float hx, float hy) { 
      return 0.36f / (hx*hy);
    }

    CUDA_CALLABLE static float kernel(float arg) {
      arg = fabsf(arg);
      float cubic = fmaxf(2.0f-arg, 0.0f);
      return 0.03125f*cubic*cubic*cubic*(arg*arg+6.0f*arg+4.0f);
    }
};

template<> 
class EpanechnikovProductConv2dTraits<double> {

  public:
    CUDA_CALLABLE static double normalization(double hx, double hy) { 
      return 0.36 / (hx*hy);
    }

    CUDA_CALLABLE static double kernel(double arg) {
      arg = fabs(arg);
      double cubic = fmax(2.0-arg, 0.0);
      return 0.03125*cubic*cubic*cubic*(arg*arg+6.0*arg+4.0);
    }
};


// GaussianConvKernel1d
// ------------------------

template<typename T> 
class GaussianConv1dTraits {

  public:
    // evaluates 1 / (sqrt(4pi)) / hx; 1 / sqrt(4pi) = 0.2820947917738781434740397257 
    static T normalization(T hx) { 
      return 0.2820947917738781434 / hx;
    }

    // evaluates exp(-0.25 * point_arg)
    static T kernel(T point_arg) {
      return std::exp(-0.25*point_arg);
    }
};

template<> 
class GaussianConv1dTraits<float> {

  public:
    CUDA_CALLABLE static float normalization(float hx) { 
      return 0.2820947917f / hx;
    }

    CUDA_CALLABLE static float kernel(float point_arg) {
      return expf(-0.25f*point_arg);
    }
};

template<> 
class GaussianConv1dTraits<double> {

  public:
    CUDA_CALLABLE static double normalization(double hx) { 
      return 0.2820947917738781434 / hx;
    }

    CUDA_CALLABLE static double kernel(double point_arg) {
      return exp(-0.25*point_arg);
    }
};

// GaussianProductConv2d
// ---------------------

template<typename T> 
class GaussianProductConv2dTraits {

  public:
    // evaluates 1 / (4pi*hx*hy); 1 / (4pi) = 0.079577471545947667884441881686
    static T normalization(T hx, T hy) { 
      return 0.0795774715459476678 / (hx*hy);
    }

    // evaluates exp(-0.5 * point_arg)
    static T kernel(T point_arg) {
      return exp(-0.25 * point_arg);
    }
};

template<> 
class GaussianProductConv2dTraits<float> {

  public:
    CUDA_CALLABLE static float normalization(float hx, float hy) { 
      return 0.0795774715f / (hx*hy);
    }

    CUDA_CALLABLE static float kernel(float point_arg) {
      return expf(-0.5f * point_arg);
    }
};

template<> 
class GaussianProductConv2dTraits<double> {

  public:
    CUDA_CALLABLE static double normalization(double hx, double hy) { 
      return 0.0795774715459476678 / (hx*hy);
    }

    CUDA_CALLABLE static double kernel(double point_arg) {
      return exp(-0.5 * point_arg);
    }
};

}

#endif
