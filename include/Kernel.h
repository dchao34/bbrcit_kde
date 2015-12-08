#ifndef BBRCITKDE_KERNEL_H__
#define BBRCITKDE_KERNEL_H__

#include <cmath>
#include <memory>

#include <Constants.h>
#include <Point.h>

// Abstract base class for Kernels. 
class Kernel1d {
  public:
    Kernel1d() {}
    virtual ~Kernel1d() {};

    std::shared_ptr<Kernel1d> clone() const { return cloneImpl(); }

    virtual double operator()(const Point1d&) = 0;

  private:
    virtual std::shared_ptr<Kernel1d> cloneImpl() const = 0;
};

// Gaussian kernel in 1 dimension: K(x) = 1/sqrt(2pi) * exp(-0.5 x*x)
class GaussKernel1d : public Kernel1d {
  public:
    GaussKernel1d() {}
    ~GaussKernel1d() {}

    inline double operator()(const Point1d&) override;
    
  private:
    std::shared_ptr<Kernel1d> cloneImpl() const override { 
      return std::shared_ptr<GaussKernel1d>(new GaussKernel1d(*this));
    }
};

double GaussKernel1d::operator()(const Point1d &p) {
  return 1.0 / kSqrt2pi * std::exp(-0.5*p[0]*p[0]);
}

#endif
