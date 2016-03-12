#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <utility>

#include <Kernels/EpanechnikovKernel.h>
#include <Kernels/EpanechnikovProductKernel2d.h>
#include <Kernels/GaussianKernel.h>
#include <KernelDensity.h>

#include "kde_test_utils.h"

using namespace std;

namespace {
  using FloatType = double;
  using KFloatType = float;
  using KernelType = bbrcit::EpanechnikovProductKernel2d<KFloatType>;
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
  int n_references = 100000;
  cout << "+ generating " << n_references << " reference points " << endl;

  default_random_engine e;
  vector<DataPointType> references;

  // weighted 
  start = std::chrono::high_resolution_clock::now();
  generate_bimodal_gaussian(e, references, n_references, 
                            1, 1, 0.5, 0.3, 30, 
                            -1, -1, 0.5, 0.3, -30, 
                            0.25, 3);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;

  fout.open("test_kde17_data.csv");
  write_scatter_data(fout, references);
  fout.close();

  cout << endl;

  // 2. build the kde's
  cout << "+ building kde's " << endl;

  size_t leaf_max = 1024;

  start = std::chrono::high_resolution_clock::now();
  KernelDensityType kde(references, leaf_max); 
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;

  cout << endl;
  
  // 3. cross validate

  cout << "+ likelihood cross validation: " << endl;

  vector<pair<double, double>> bandwidths_pair;
  double start_x = 0.38, end_x = 0.45; int steps_x = 4;
  double start_y = 0.11, end_y = 0.15; int steps_y = 4;
  double delta_x = (end_x-start_x)/steps_x;
  double delta_y = (end_y-start_y)/steps_y;
  for (int j = 0; j < steps_y; ++j) {
    for (int i = 0; i < steps_x; ++i) {
      bandwidths_pair.push_back({start_x+i*delta_x, start_y+j*delta_y});
    }
  }
  
  FloatType cv, best_cv = -std::numeric_limits<FloatType>::max();
  double best_bw_x, best_bw_y;
  for (const auto &p : bandwidths_pair) {
    kde.kernel().set_bandwidths(p.first, p.second);
    cv = kde.likelihood_cross_validate();
    if (cv > best_cv) { 
      best_cv = cv; 
      best_bw_x = p.first; best_bw_y = p.second; 
    }
    std::cout << "  (" << p.first << ", " << p.second << ") " << cv << std::endl;
  }

  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  cout << std::endl;

  cout << "  best bandwidth: (" << best_bw_x << ", " << best_bw_y << "), cv score = " << best_cv << std::endl;

  cout << std::endl;

#ifndef __CUDACC__
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;
#else 
  cout << "  gpu time: " << elapsed.count() << " ms. " << std::endl;
#endif

  cout << std::endl;
  
 
  // 4. generate the query grid
  vector<DataPointType> grid, queries;
  start_x = -2, end_x = 2; steps_x = 100;
  start_y = -2, end_y = 2; steps_y = 100;

  cout << "+ generating " << steps_x << "x" << steps_y << " query grid" << endl;

  start = std::chrono::high_resolution_clock::now();
  generate_2dgrid(grid, start_x, end_x, steps_x, start_y, end_y, steps_y);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;

  cout << endl;

  // 5. evaluate cross validated density
  cout << "+ evaluating cross validated density" << endl; 
  queries = grid;

  kde.kernel().set_bandwidths(0.38, 0.12);

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

  fout.open("test_kde17_cv.csv");
  write_kde2d_result(fout, queries, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);
  fout.close();

  cout << endl;
  
  // 6. adaptive density
  cout << "+ adapting density" << endl; 
  start = std::chrono::high_resolution_clock::now();
  kde.adapt_density(0.5);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
#ifndef __CUDACC__
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;
#else 
  cout << "  gpu time: " << elapsed.count() << " ms. " << std::endl;
#endif

  cout << endl;

  // 7. evaluate adapted density
  cout << "+ evaluating adapted density" << endl; 
  queries = grid;

  start = std::chrono::high_resolution_clock::now();
  kde.eval(queries, rel_tol, abs_tol, leaf_max);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
#ifndef __CUDACC__
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;
#else 
  cout << "  gpu time: " << elapsed.count() << " ms. " << std::endl;
#endif

  fout.open("test_kde17_adapt.csv");
  write_kde2d_result(fout, queries, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);
  fout.close();

  cout << endl;
  


  return 0;
}
