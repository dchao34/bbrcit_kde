#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <thread>
#include <mutex>

#include <KernelDensity.h>
#include <Kernels/GaussianKernel.h>
#include <Kernels/EpanechnikovKernel.h>
#include <Point.h>

#include "kde_test_utils.h"

using namespace std;

using bbrcit::KernelDensity;
using bbrcit::EpanechnikovKernel;
using bbrcit::GaussianKernel;
using bbrcit::Point;

using EpanKdeType = KernelDensity<2,EpanechnikovKernel<2>>;
using DataPointType = typename EpanKdeType::DataPointType;

mutex m;

void single_tree_evaluate_segment(
    int b, int e, vector<DataPointType> &queries, 
    double rel_err, double abs_err, 
    const EpanKdeType *epan_kde) {
  
  for (int i = b; i < e; ++i) {
    auto result = epan_kde->eval(queries[i], rel_err, abs_err); 
    auto attr = queries[i].attributes();
    attr.set_lower(result);
    attr.set_upper(result);
    queries[i].set_attributes(attr);
  }
}

int main() {

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed;

  mt19937 e;
  normal_distribution<> gaussian(0, 1);

  vector<DataPointType> data, queries;
  double rel_err = 1e-6; double abs_err = 1e-10;

  int k = 17;
  int n_samples = std::pow(2, k);

  for (int i = 0; i < n_samples; ++i) { 
    data.push_back({{gaussian(e), gaussian(e)}}); 
  }

  // kde build
  start = std::chrono::high_resolution_clock::now();
  EpanKdeType *epan_kde = new EpanKdeType(data); 
  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "kde build time: " << elapsed.count() << " " << endl;;

  // set bandwidth
  double epan_bw = 0.1;
  epan_kde->kernel().set_bandwidth(epan_bw);

  // single tree evaluation
  queries = data;
  start = std::chrono::high_resolution_clock::now();

  single_tree_evaluate_segment(0, queries.size(), queries, 
                               rel_err, abs_err, epan_kde);

  end = std::chrono::high_resolution_clock::now();

  elapsed = end-start;
  cout << "single tree time: " << elapsed.count() << endl;

  // single tree parallel evaluation
  queries = data;
  start = std::chrono::high_resolution_clock::now();

  int n_workers = 2;
  //int n_workers = thread::hardware_concurrency();
  vector<thread> children(n_workers-1);

  int workload = queries.size() / n_workers;
  int first = 0;
  for (int i = 0; i < n_workers-1; ++i) {
    children[i] = thread(single_tree_evaluate_segment, 
                         first, first + workload, 
                         ref(queries), rel_err, abs_err, epan_kde); 
    first += workload;
  }
  single_tree_evaluate_segment(first, queries.size(), queries, 
                               rel_err, abs_err, epan_kde);

  for (auto &t : children) { t.join(); }

  end = std::chrono::high_resolution_clock::now();

  elapsed = end-start;
  cout << "parallel single tree time: " << elapsed.count() << endl;

  delete epan_kde;


  return 0;
}
