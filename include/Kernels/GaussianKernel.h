#ifndef BBRCITKDE_GAUSSIANKERNEL_H__
#define  BBRCITKDE_GAUSSIANKERNEL_H__

#define _USE_MATH_DEFINES

#ifdef __CUDACC__
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_CALLABLE
#endif

#include <cmath>
#include <iostream>

#include <KdeTraits.h>
#include <Point2d.h>

namespace bbrcit {

// The Gaussian kernel in D dimensions of bandwidth h is
// defined as follows:
//
// K(x) = unit_volume / h^D * exp{ -0.5 * x'x / h*h}
//
// where: 
//
// + x'x: 
//     dot product of x. 
// + unit_volume: 
//     (2 \pi)^{-D/2}
//
// The two point version is defined as follows:
//
// K(x, y) := K(x-y)
//
template<int D, typename T=double>
class GaussianKernel {

  public:

    static constexpr int dim() { return D; }

  public:
    CUDA_CALLABLE GaussianKernel();
    CUDA_CALLABLE GaussianKernel(T bandwidth);

    CUDA_CALLABLE ~GaussianKernel() = default;
    CUDA_CALLABLE GaussianKernel(const GaussianKernel&) = default;
    CUDA_CALLABLE GaussianKernel(GaussianKernel&&) = default;
    CUDA_CALLABLE GaussianKernel& operator=(const GaussianKernel&) = default;
    CUDA_CALLABLE GaussianKernel& operator=(GaussianKernel&&) = default;

    // compute the normalization
    CUDA_CALLABLE T normalization() const;

    // evaluate the two point kernel, but do not include the 
    // normalization factor. e.g. evaluates exp{-0.5 * (x-y)'(x-y)/(h*h)
    template<typename PointT>
    CUDA_CALLABLE T unnormalized_eval(const PointT&, const PointT&) const;

    // get/set the bandwidth
    CUDA_CALLABLE T bandwidth() const;
    CUDA_CALLABLE void set_bandwidth(T);

  private:
    T bandwidth_;

    template<typename PointT>
    CUDA_CALLABLE T point_arg_eval(const PointT&, const PointT&) const;
    
    // TODO: allow these overloads for now. would be better if only some are 
    // callable when D has the right dimension. 
    // consider using IfThenElse metafunction to custimize these overloads?
    CUDA_CALLABLE T point_arg_eval(const Point2d<T>&, const Point2d<T>&) const;
    
};

// Implementations
// ---------------

template<int D, typename T>
GaussianKernel<D, T>::GaussianKernel() : bandwidth_(1.0) {}

template<int D, typename T>
GaussianKernel<D,T>::GaussianKernel(T bw) : bandwidth_(bw) {}

template<int D, typename T>
inline T GaussianKernel<D,T>::bandwidth() const { return bandwidth_; }

template<int D, typename T>
inline void GaussianKernel<D,T>::set_bandwidth(T bw) { bandwidth_ = bw; }

template<int D, typename T>
  template<typename PointT>
T GaussianKernel<D,T>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return exp(-0.5 * point_arg_eval(p, q));
}

template <int D, typename T>
inline T GaussianKernel<D,T>::normalization() const {
  return pow(2*M_PI, -D/2.0) / pow(bandwidth_, D);
}

// evaluates the (x-y)'(x-y)/h*h part for the kernel argument. 
template<int D, typename T>
  template<typename PointT>
T GaussianKernel<D,T>::point_arg_eval(const PointT &lhs, const PointT &rhs) const {

  T result = ConstantTraits<T>::zero(); 

  T diff; 
  for (int i = 0; i < D; ++i) {
    diff = lhs[i] - rhs[i];
    result += diff * diff;
  }

  return result / (bandwidth_ * bandwidth_);
}

template<int D, typename T>
inline T GaussianKernel<D,T>::point_arg_eval(
    const Point2d<T> &lhs, 
    const Point2d<T> &rhs) const {
  return ((lhs.x()-rhs.x())*(lhs.x()-rhs.x()) +
          (lhs.y()-rhs.y())*(lhs.y()-rhs.y())) / (bandwidth_*bandwidth_);
}


// Specializations
// ---------------

// 1. normalization()
// ------------------

template <>
inline float GaussianKernel<1,float>::normalization() const {
  // 1 / (sqrt(2pi) * h); 1 / sqrt(2pi) = 0.39894228040143267794
  return 0.3989422804f / bandwidth_;
}

template <>
inline double GaussianKernel<1,double>::normalization() const {
  // 1 / (sqrt(2pi) * h); 1 / sqrt(2pi) = 0.39894228040143267794
  return 0.3989422804014326779 / bandwidth_;
}

template <>
inline float GaussianKernel<2,float>::normalization() const {
  // 1 / (2pi) * h * h; 1 / (2pi) = 0.15915494309189533577
  return 0.1591549431f / (bandwidth_ * bandwidth_);
}

template <>
inline double GaussianKernel<2,double>::normalization() const {
  // 1 / (2pi) * h * h; 1 / (2pi) = 0.15915494309189533577
  return 0.1591549430918953357 / (bandwidth_ * bandwidth_);
}

// 2. unnormalized_eval()
// ----------------------

template<>
  template<typename PointT>
inline float GaussianKernel<1,float>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return expf(-0.5f * point_arg_eval(p, q));
}

template<>
  template<typename PointT>
inline double GaussianKernel<1,double>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return exp(-0.5 * point_arg_eval(p, q));
}

template<>
  template<typename PointT>
inline float GaussianKernel<2,float>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return expf(-0.5f * point_arg_eval(p, q));
}

template<>
  template<typename PointT>
inline double GaussianKernel<2,double>::unnormalized_eval(const PointT &p, const PointT &q) const {
  return exp(-0.5 * point_arg_eval(p, q));
}

}

#endif
