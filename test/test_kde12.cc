#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>

#include <Kernels/EpanechnikovProductKernel2d.h>
#include <Kernels/GaussianProductKernel2d.h>
#include <KernelDensity.h>

#include "kde_test_utils.h"

using namespace std;

namespace {
  using FloatType = double;
  using KFloatType = float;
  using KernelType = bbrcit::EpanechnikovProductKernel2d<KFloatType>;
  //using KernelType = bbrcit::GaussianProductKernel2d<KFloatType>;
  using KernelDensityType = bbrcit::KernelDensity<2, KernelType, FloatType>;
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
  generate_bimodal_gaussian(e, references, n_references, 
                            1, 1, 0.5, 0.3, 30, 
                            -1, -1, 0.5, 0.3, -30);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;

  fout.open("test_kde12_data.csv");
  write_scatter_data(fout, references);
  fout.close();

  cout << endl;

  // 2. generate the query grid
  vector<DataPointType> grid, queries;
  double start_x = -2, end_x = 2; int steps_x = 100;
  double start_y = -2, end_y = 2; int steps_y = 100;

  cout << "+ generating " << steps_x << "x" << steps_y << " query grid" << endl;

  start = std::chrono::high_resolution_clock::now();
  generate_2dgrid(grid, start_x, end_x, steps_x, start_y, end_y, steps_y);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;

  cout << endl;

  // 3. build the kernel density estimator
  cout << "+ building kde (kdtree construction)" << endl;

  size_t leaf_max = 1024;

  start = std::chrono::high_resolution_clock::now();
  KernelDensityType kde(references, leaf_max); 
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;
  
  // configure the kernel
  kde.kernel().set_bandwidths(0.2, 0.05);

  cout << endl;
  

  // 4. direct kde evaluation
  cout << "+ direct evaluation" << endl; 
  queries = grid;

  start = std::chrono::high_resolution_clock::now();
  kde.direct_eval(queries);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
#ifndef __CUDACC__
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;
#else 
  cout << "  gpu time: " << elapsed.count() << " ms. " << std::endl;
#endif

#ifndef __CUDACC__
  fout.open("test_kde12_cpu_direct.csv");
#else 
  fout.open("test_kde12_gpu_direct.csv");
#endif

  write_kde2d_result(fout, queries, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);
  fout.close();

  cout << endl;

  // 4. dual tree evaluation
  cout << "+ dual tree evaluation" << endl; 
  queries = grid;

  FloatType rel_tol = 1e-6, abs_tol = 1e-6;

  start = std::chrono::high_resolution_clock::now();
  kde.eval(queries, rel_tol, abs_tol, leaf_max);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
#ifndef __CUDACC__
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;
#else 
  cout << "  gpu time: " << elapsed.count() << " ms. " << std::endl;
#endif

#ifndef __CUDACC__
  fout.open("test_kde12_cpu_tree.csv");
#else 
  fout.open("test_kde12_gpu_tree.csv");
#endif
  write_kde2d_result(fout, queries, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);
  fout.close();

  cout << endl;


  return 0;
}
