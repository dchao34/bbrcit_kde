#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>

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

int main() {

  using KernelDensityType = KernelDensity<1,GaussianKernel<1>>;
  //using KernelDensityType = KernelDensity<1,EpanechnikovKernel<1>>;
  using DataPointType = typename KernelDensityType::DataPointType;

  mt19937 e;
  normal_distribution<> gaussian(0, 1);

  vector<DataPointType> data, queries;
  double rel_err = 1e-6; double abs_err = 1e-12;

  int n_samples = 10000;
  for (int i = 0; i < n_samples; ++i) { data.push_back({{gaussian(e)}}); }

  int n_steps = 100;
  generate_1dgrid(queries, -5, 5, n_steps);

  KernelDensityType kde(data); 
  kde.kernel().set_bandwidth(0.1);

  for (auto &p : queries) { 
    auto result = kde.eval(p, rel_err, abs_err); 
    auto attr = p.attributes();
    attr.set_lower(result);
    attr.set_upper(result);
    p.set_attributes(attr);
  }


  return 0;
}
