#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>

#include <Kernels/EpanechnikovProductKernel2d.h>
#include <Kernels/EpanechnikovKernel.h>
#include <KernelDensity.h>

#include <MarginalDensityMaker.h>

#include "kde_test_utils.h"

using namespace std;

namespace {
  using FloatType = double;
  using KFloatType = float;

  double rel_tol = 1e-6;
  double abs_tol = 1e-8;

  using Kernel2dType = bbrcit::EpanechnikovProductKernel2d<KFloatType>;
  using KernelDensity2dType = bbrcit::KernelDensity<2, Kernel2dType, FloatType>;
  using DataPoint2dType = typename KernelDensity2dType::DataPointType;

  using Kernel1dType = bbrcit::EpanechnikovKernel<1, KFloatType>;
  using KernelDensity1dType = bbrcit::KernelDensity<1, Kernel1dType, FloatType>;
  using DataPoint1dType = typename KernelDensity1dType::DataPointType;
}

int main() {

#ifdef __CUDACC__
  cudaSetDevice(1);
  cudaDeviceSynchronize();
#endif

  ofstream fout;

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double, std::milli> elapsed;

  // 1. generate the reference points
  int n_references = 10000;
  cout << "+ generating " << n_references << " reference points " << endl;

  default_random_engine e;
  vector<DataPoint2dType> references;

  start = std::chrono::high_resolution_clock::now();
  generate_bimodal_gaussian(e, references, n_references, 
                            1, 1, 0.5, 0.3, 30, 
                            -1, -1, 0.5, 0.3, -30);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  running time: " << elapsed.count() << " ms. " << std::endl;

  fout.open("test1d_data.csv");
  write_scatter_data(fout, references);
  fout.close();

  cout << endl;

  // 2. generate the query grid
  vector<DataPoint2dType> grid2d;
  double start_x = -3, end_x = 2; int steps_x = 100;
  double start_y = -2, end_y = 2.5; int steps_y = 100;

  cout << "+ generating " << steps_x << "x" << steps_y << " query grid" << endl;

  start = std::chrono::high_resolution_clock::now();
  generate_2dgrid(grid2d, start_x, end_x, steps_x, start_y, end_y, steps_y);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  running time: " << elapsed.count() << " ms. " << std::endl;

  cout << endl;

  // 3. build the kernel density estimator
  cout << "+ building kde (kdtree construction)" << endl;

  size_t leaf_max = 1024;

  start = std::chrono::high_resolution_clock::now();
  KernelDensity2dType kde(references, leaf_max); 
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  running time: " << elapsed.count() << " ms. " << std::endl;
  
  // configure the kernel
  kde.kernel().set_bandwidths(0.01, 0.2);

  cout << endl;

  // 4. dual tree evaluation
  cout << "+ dual tree evaluation" << endl; 

  start = std::chrono::high_resolution_clock::now();
  kde.eval(grid2d, rel_tol, abs_tol, leaf_max);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  running time: " << elapsed.count() << " ms. " << std::endl;

  fout.open("test1d.csv");
  write_kde2d_result(fout, grid2d, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);
  fout.close();

  cout << endl;


  // marginal in x
  std::vector<DataPoint1dType> grid1d;
  generate_1dgrid(grid1d, start_x, end_x, steps_x);

  auto kde_x = bbrcit::EpanechnikovProduct2dKdeMarginalizer<KFloatType,FloatType>::make_marginal_density_x(kde);
  kde_x.eval(grid1d, rel_tol, abs_tol, leaf_max);
  fout.open("test1d_x.csv");
  write_kde1d_result(fout, grid1d);
  fout.close();

  generate_1dgrid(grid1d, start_y, end_y, steps_y);
  auto kde_y = bbrcit::EpanechnikovProduct2dKdeMarginalizer<KFloatType,FloatType>::make_marginal_density_y(kde);
  kde_y.eval(grid1d, rel_tol, abs_tol, leaf_max);
  fout.open("test1d_y.csv");
  write_kde1d_result(fout, grid1d);
  fout.close();

  return 0;
}
