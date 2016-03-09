#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>

#include <Kernels/GaussianKernel.h>
#include <Kernels/EpanechnikovKernel.h>
#include <KernelDensity.h>
#include <Point.h>

#include "kde_test_utils.h"

using namespace std;

namespace {
  using FloatType = double;
  using KFloatType = float;
  //using KernelType = bbrcit::EpanechnikovKernel<2, KFloatType>;
  using KernelType = bbrcit::GaussianKernel<2, KFloatType>;
  using KernelDensityType = bbrcit::KernelDensity<2, KernelType, FloatType>;
  using DataPointType = typename KernelDensityType::DataPointType;
}

void print_line(
    std::ostream &os, 
    int k, double direct_time, double tree_time) {
  os << k << " ";
  os << direct_time << " ";
  os << tree_time << " ";
  os << std::endl;
}

int main() {

#ifdef __CUDACC__
  cudaDeviceSynchronize();
#endif

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double, std::milli> elapsed;

  mt19937 e;
  normal_distribution<FloatType> d(0, 1);

  vector<DataPointType> data, queries;
  double rel_err = 1e-6, abs_err = 1e-10;

#ifndef __CUDACC__
  ofstream fout("test_kde11_cpu.csv");
#else 
  ofstream fout("test_kde11_gpu.csv");
#endif

#ifndef __CUDACC__
  int leaf_max = 32;
#else 
  int leaf_max = 4096;
  int block_size = 128;
#endif

#ifndef __CUDACC__
  int kstart = 7;
  int kend = 20;
  int direct_kmax = 13;
#else 
  int kstart = 7;
  int kend = 21;
  int direct_kmax = 18;
#endif

  size_t n_samples; double direct_time, tree_time;
  for (int k = kstart; k <= kend; ++k) {

    data.clear();
    n_samples = 2 << k;
    for (size_t i = 0; i < n_samples; ++i) {
      data.push_back({{d(e), d(e)}});
    }

    KernelDensityType *kde = new KernelDensityType(data, leaf_max);
    kde->kernel().set_bandwidth(0.1);

    if (k <= direct_kmax) {
      queries = data;

      start = std::chrono::high_resolution_clock::now();
      kde->direct_eval(queries);
      end = std::chrono::high_resolution_clock::now();
      elapsed = end - start;

      direct_time = elapsed.count();

    } else {

      direct_time = 0;

    }

    queries = data;

    start = std::chrono::high_resolution_clock::now();
#ifndef __CUDACC__
    kde->eval(queries, rel_err, abs_err, leaf_max);
#else 
    kde->eval(queries, rel_err, abs_err, leaf_max, block_size);
#endif
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;

    tree_time = elapsed.count();

    print_line(cout, k, direct_time, tree_time);
    print_line(fout, k, direct_time, tree_time);

    delete kde;


  }

  return 0;
}
