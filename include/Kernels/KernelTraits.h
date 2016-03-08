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

    // evaluates 1 / (2pi) * h * h; 1 / (2pi) = 0.15915494309189533577
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

    // evaluates 1 / (2pi) * h * h; 1 / (2pi) = 0.15915494309189533577
    CUDA_CALLABLE static double normalization(double h) { 
      return 0.1591549430918953357 / (h*h);
    }

    CUDA_CALLABLE static double kernel(double point_arg) {
      return exp(-0.5 * point_arg);
    }
};

}

#endif
