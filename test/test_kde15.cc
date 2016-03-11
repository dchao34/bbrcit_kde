#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <limits>

#include <Kernels/EpanechnikovKernel.h>
#include <Kernels/GaussianKernel.h>
#include <KernelDensity.h>

#include "kde_test_utils.h"

using namespace std;

namespace {
  using FloatType = double;
  using KFloatType = float;
  using KernelType = bbrcit::EpanechnikovKernel<1, KFloatType>;
  //using KernelType = bbrcit::GaussianKernel<1, KFloatType>;
  using KernelDensityType = bbrcit::KernelDensity<1, KernelType, FloatType>;
  using DataPointType = typename KernelDensityType::DataPointType;
}

int main() {

#ifdef __CUDACC__
  cudaDeviceSynchronize();
#endif

  ofstream fout;

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double, std::milli> elapsed;

  // 1. generate the reference points
  int n_references = 10000;
  cout << "+ generating " << n_references << " reference points " << endl;

  default_random_engine e;
  vector<DataPointType> references;

  start = std::chrono::high_resolution_clock::now();
  generate_bimodal_gaussian(e, references, n_references, 0.3, 0.1, -0.3, 0.1, 0.75);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;

  cout << endl;

  // 2. build the kernel density estimator
  cout << "+ building kde (kdtree construction)" << endl;

  size_t leaf_max = 1024;

  start = std::chrono::high_resolution_clock::now();
  KernelDensityType kde(references, leaf_max); 
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;
  
  // configure the kernel
  kde.kernel().set_bandwidth(0.2);

  cout << endl;
  
  // 3. cross validate
  cout << "+ likelihood cross validation. " << endl;
  cout << std::endl;

  start = std::chrono::high_resolution_clock::now();

  std::vector<FloatType> bandwidths = { 0.0478, 0.048, 0.049, 0.05, 0.06, 0.07, 0.1, 0.2 };

  FloatType best_bw;
  FloatType cv, best_cv = std::numeric_limits<FloatType>::min();
  for (size_t i = 0; i < bandwidths.size(); ++i) {
    kde.kernel().set_bandwidth(bandwidths[i]);
    cv = kde.likelihood_cross_validate();
    if (cv > best_cv) { best_cv = cv; best_bw = bandwidths[i]; }
    std::cout << "  " << kde.kernel().bandwidth() << " " << cv << std::endl;
  }

  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  cout << std::endl;

  cout << "  best bandwidth: " << best_bw << ", cv score = " << best_cv << std::endl;

  cout << std::endl;

#ifndef __CUDACC__
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;
#else 
  cout << "  gpu time: " << elapsed.count() << " ms. " << std::endl;
#endif

  return 0;
}
