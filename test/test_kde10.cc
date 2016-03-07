#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>

#include <Kernels/EpanechnikovKernel.h>
#include <KernelDensity.h>
#include <Point.h>

#include "kde_test_utils.h"

using namespace std;

namespace {
  using FloatType = double;
  using KernelType = bbrcit::EpanechnikovKernel<2, FloatType>;
  using KernelDensityType = bbrcit::KernelDensity<2, FloatType, KernelType>;
  using DataPointType = typename KernelDensityType::DataPointType;
}

int main() {

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

  fout.open("test_kde10_data.csv");
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

  start = std::chrono::high_resolution_clock::now();
  KernelDensityType kde(references, 2); 
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;
  
  // configure the kernel
  kde.kernel().set_bandwidth(0.2);

  cout << endl;
  

  // 4. direct kde evaluation
  cout << "+ direct kde evaluation" << endl; 
  queries = grid;

#ifdef __CUDACC__
  cudaDeviceSynchronize();
#endif

  start = std::chrono::high_resolution_clock::now();
  kde.direct_eval(queries);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cout << "  cpu time: " << elapsed.count() << " ms. " << std::endl;

  fout.open("test_kde10_direct.csv");
  write_kde2d_result(fout, queries, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);
  fout.close();

  cout << endl;

  return 0;
}
