#include <iostream>
#include <stdio.h>

#include <Point2d.h>
#include <Kernels/GaussianKernel.h>

#ifdef __CUDACC__
__global__ void kernel() {
  Point2d<float> p(1,2,3);
  printf("%f, %f, %f\n", p.x(), p.y(), p.w());
  return;
}
#endif

int main() {

  Point2d<float> p;
  std::cout << sizeof(p) << std::endl;

#ifdef __CUDACC__
  kernel<<<1,1>>>();
  cudaDeviceSynchronize();
#endif

  return 0;
}
