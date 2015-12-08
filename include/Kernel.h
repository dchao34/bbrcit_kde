#ifndef BBRCITKDE_KERNEL_H__
#define BBRCITKDE_KERNEL_H__

#include <cmath>

#include <Constants.h>
#include <Point.h>

class Kernel {
  public:
    Kernel() {}
    virtual ~Kernel() = 0;
};

Kernel::~Kernel() {}

class GaussKernel1d : public Kernel {
  public:
    GaussKernel1d() {}
    ~GaussKernel1d() {}
    double inline operator()(const Point1d&);
};

double GaussKernel1d::operator()(const Point1d &p) {
  return 1.0 / kSqrt2pi * std::exp(-0.5*p[0]*p[0]);
}

#endif
