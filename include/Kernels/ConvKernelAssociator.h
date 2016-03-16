#ifndef BBRCITKDE_CONVKERNELASSOCIATOR_H__
#define BBRCITKDE_CONVKERNELASSOCIATOR_H__

#include <Kernels/EpanechnikovKernel.h>
#include <Kernels/EpanechnikovConvKernel1d.h>
#include <Kernels/EpanechnikovProductKernel2d.h>
#include <Kernels/EpanechnikovProductConvKernel2d.h>
#include <Kernels/GaussianKernel.h>
#include <Kernels/GaussianConvKernel1d.h>
#include <Kernels/GaussianProductKernel2d.h>
#include <Kernels/GaussianProductConvKernel2d.h>

namespace bbrcit {

// class template that offers meta functions to associate kernels 
// with its convolution kernel. 
template<typename T> 
class ConvKernelAssociator {};

// Epanechnikov Kernels
// --------------------

template<typename T>
class ConvKernelAssociator<EpanechnikovKernel<1,T>> {

  public:

    // kernel type
    using KernelType = EpanechnikovKernel<1,T>;

    // convolution kernel type associated with the kernel
    using ConvKernelType = EpanechnikovConvKernel1d<T>;

    // given an object of `KernelType`, return a (copy of) `ConvKernelType` object.
    static ConvKernelType make_convolution_kernel(const KernelType &k) {
      return ConvKernelType(k.bandwidth());
    }
};

template<typename T>
class ConvKernelAssociator<EpanechnikovProductKernel2d<T>> {

  public:

    using KernelType = EpanechnikovProductKernel2d<T>;
    using ConvKernelType = EpanechnikovProductConvKernel2d<T>;

    static ConvKernelType make_convolution_kernel(const KernelType &k) {
      return ConvKernelType(k.hx(), k.hy());
    }
};

// Gaussian Kernels
// ----------------

template<typename T>
class ConvKernelAssociator<GaussianKernel<1,T>> {

  public:

    // kernel type
    using KernelType = GaussianKernel<1,T>;

    // convolution kernel type associated with the kernel
    using ConvKernelType = GaussianConvKernel1d<T>;

    // given an object of `KernelType`, return a (copy of) `ConvKernelType` object.
    static ConvKernelType make_convolution_kernel(const KernelType &k) {
      return ConvKernelType(k.bandwidth());
    }
};

template<typename T>
class ConvKernelAssociator<GaussianProductKernel2d<T>> {

  public:

    using KernelType = GaussianProductKernel2d<T>;
    using ConvKernelType = GaussianProductConvKernel2d<T>;

    static ConvKernelType make_convolution_kernel(const KernelType &k) {
      return ConvKernelType(k.hx(), k.hy());
    }
};

}

#endif
