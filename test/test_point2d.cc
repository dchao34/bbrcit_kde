#include <iostream>
#include <iomanip>
#include <stdio.h>

#include <Point2d.h>
#include <Kernels/EpanechnikovKernel.h>
#include <Kernels/GaussianKernel.h>

using FloatType = double;

using bbrcit::Point2d;
using bbrcit::GaussianKernel;
using bbrcit::EpanechnikovKernel;

#ifdef __CUDACC__
__global__ void kernel() {
  Point2d<FloatType> origin, p(1,1,1);
  GaussianKernel<2,FloatType> k(1.0);
  printf("%1.12f (c.f. 0.05854983152431917)\n", k.normalization() * k.unnormalized_eval(p, origin));
  k.set_bandwidth(2.0);
  printf("%1.12f (c.f. 0.03098749857741324)\n", k.normalization() * k.unnormalized_eval(p, origin));

  EpanechnikovKernel<2,FloatType> ek(1.0);
  Point2d<FloatType> q(3,2,1);
  printf("%1.12f (c.f. 0.0)\n", ek.normalization() * ek.unnormalized_eval(p, q));
  ek.set_bandwidth(5.0);
  printf("%1.12f (c.f. 0.020371832715762605)\n", ek.normalization() * ek.unnormalized_eval(p, q));

  return;
}
#endif

int main() {

  std::cout << "sizeof(Point2d<>) = " << sizeof(Point2d<FloatType>) << std::endl;
  std::cout << "sizeof(EpanechnikovKernel<>) = " << sizeof(EpanechnikovKernel<2,FloatType>) << std::endl;
  std::cout << "sizeof(GaussianKernel<>) = " << sizeof(GaussianKernel<2,FloatType>) << std::endl;
  std::cout << std::endl;

  std::cout << std::setprecision(12) << std::fixed;

  std::cout << "cpu tests: " << std::endl;
  Point2d<FloatType> origin, p(1,1,1);

  GaussianKernel<2,FloatType> k(1.0);
  std::cout << (k.normalization() * k.unnormalized_eval(p, origin));
  std::cout << " (c.f. 0.05854983152431917)" << std::endl;
  k.set_bandwidth(2.0);
  std::cout << (k.normalization() * k.unnormalized_eval(p, origin));
  std::cout << " (c.f. 0.03098749857741324)" << std::endl;

  EpanechnikovKernel<2,FloatType> ek(1.0);
  Point2d<FloatType> q(3,2,1);
  std::cout << (ek.normalization() * ek.unnormalized_eval(p, q));
  std::cout << " (c.f. 0.0)" << std::endl;
  ek.set_bandwidth(5.0);
  std::cout << (ek.normalization() * ek.unnormalized_eval(p, q));
  std::cout << " (c.f. 0.020371832715762605)" << std::endl;
  std::cout << std::endl;

#ifdef __CUDACC__
  std::cout << "gpu tests: " << std::endl;
  kernel<<<1,1>>>();
  cudaDeviceSynchronize();
  std::cout << std::endl;
#endif

  return 0;
}
