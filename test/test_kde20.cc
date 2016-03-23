#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include <random>
#include <cmath>

#include <Kernels/EpanechnikovProductKernel2d.h>
#include <KernelDensity.h>

namespace {
  using KernelFloatType = float;
  using FloatType = double;
  using KernelType = bbrcit::EpanechnikovProductKernel2d<KernelFloatType>;
  using KernelDensityType = bbrcit::KernelDensity<2,KernelType,FloatType>;
  using DataPointType = KernelDensityType::DataPointType;
}

int main() {

#ifdef __CUDACC__
  cudaSetDevice(0);
#endif

  // general parameters
  FloatType rel_tol = 1e-6, abs_tol = 1e-8;
  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double, std::milli> elapsed;
  std::mt19937 e;


  // generate reference points from standard normal
  size_t n_ref = 10000;

  std::cout << "+ generating " << n_ref << " points from standard normal. \n" << std::endl;

  std::normal_distribution<FloatType> d(0, 1);
  std::vector<DataPointType> ref_pts;
  for (size_t i = 0; i < n_ref; ++i) {
    ref_pts.push_back({{d(e), d(e)}, {1.0/n_ref}});
  }

  // configure kernel density
  size_t ref_leaf_max = 1024;

  std::cout << "+ creating kernel density. " << std::endl;

  start = std::chrono::high_resolution_clock::now();
  KernelDensityType kde(ref_pts, ref_leaf_max);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "  run time: " << elapsed.count() << " ms. \n" << std::endl;;

  kde.kernel().set_bandwidths(0.5, 0.5);

  // compute least square score with convolution kernel. 
  std::cout << "+ computing least squarse cv score using convolution kernel. " << std::endl;

  start = std::chrono::high_resolution_clock::now();

#ifndef __CUDACC__
  FloatType conv_cv = kde.leastsquares_cross_validate(rel_tol, abs_tol);
#else
  size_t block_size = 128;
  FloatType conv_cv = kde.leastsquares_cross_validate(rel_tol, abs_tol, block_size);
#endif

  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "  run time: " << elapsed.count() << " ms. \n" << std::endl;;

  // compute least square score with numeical integration. 
  FloatType start_x = -5, end_x = 5; int steps_x = 1000;
  FloatType start_y = -5, end_y = 5; int steps_y = 1000;

  std::cout << "+ computing least squarse cv score using numerical integration. " << std::endl;
  std::cout << "  integration grid: " << steps_x << " by " << steps_y << ", over ";
  std::cout << "[" << start_x << ", " << end_x << "]" << " x ";
  std::cout << "[" << start_y << ", " << end_y << "]" << std::endl;

  start = std::chrono::high_resolution_clock::now();

#ifndef __CUDACC__
  FloatType num_cv = leastsquares_numint_cross_validate(
      kde, start_x, end_x, steps_x, start_y, end_y, steps_y,
      rel_tol, abs_tol);
#else
  int query_leaf_max = 65536;
  FloatType num_cv = leastsquares_numint_cross_validate(
      kde, start_x, end_x, steps_x, start_y, end_y, steps_y,
      rel_tol, abs_tol, query_leaf_max, block_size);
#endif

  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "  run time: " << elapsed.count() << " ms. \n" << std::endl;;
 
  // output result
  std::cout << "+ lsq cv score for bandwidths ";
  std::cout << "(" << kde.kernel().hx() << ", " << kde.kernel().hy() << ")" << std::endl;

  std::cout << std::setprecision(10) << std::fixed;
  std::cout << "  ==> score by numerical integration: " << num_cv << std::endl;
  std::cout << "  ==> score by convolution kernel:    " << conv_cv << "\n " << std::endl;


  return 0;
}
