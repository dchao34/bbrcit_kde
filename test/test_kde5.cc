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

  using KernelDensityType = KernelDensity<1,EpanechnikovKernel<1>>;
  using DataPointType = typename KernelDensityType::DataPointType;

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed;

  mt19937 e;
  normal_distribution<> gaussian(0, 1);

  ofstream fout("test_kde5.csv");

  vector<DataPointType> data, queries;
  double rel_err = 1e-6; double abs_err = 1e-10;

  int n_samples;
  for (int k = 7; k < 22; ++k) {

    cout << k << " ";
    fout << k << " ";

    n_samples = std::pow(2, k);
    for (int i = 0; i < n_samples; ++i) { data.push_back({{gaussian(e)}}); }

    // kde build
    start = std::chrono::high_resolution_clock::now();
    KernelDensityType kde(data); 
    end = std::chrono::high_resolution_clock::now();
    elapsed = end-start;
    cout << elapsed.count() << " ";
    fout << elapsed.count() << " ";

    kde.kernel().set_bandwidth(0.001);

    // naive evaluation
    if (k < 13) {
      queries = data;
      start = std::chrono::high_resolution_clock::now();
      for (auto &p : queries) { 
        auto result = kde.naive_eval(p); 
        auto attr = p.attributes();
        attr.set_lower(result);
        attr.set_upper(result);
        p.set_attributes(attr);
      }
      end = std::chrono::high_resolution_clock::now();

      elapsed = end-start;
      cout << elapsed.count() << " ";
      fout << elapsed.count() << " ";
    } else {
      cout << 0 << " ";
      fout << 0 << " ";
    }

    // single tree evaluation
    queries = data;
    start = std::chrono::high_resolution_clock::now();
    for (auto &p : queries) { 
      auto result = kde.eval(p, rel_err, abs_err); 
      auto attr = p.attributes();
      attr.set_lower(result);
      attr.set_upper(result);
      p.set_attributes(attr);
    }
    end = std::chrono::high_resolution_clock::now();

    elapsed = end-start;
    cout << elapsed.count() << " ";
    fout << elapsed.count() << " ";

    // dual tree evaluation
    queries = data;

    start = std::chrono::high_resolution_clock::now();
    kde.eval(queries, rel_err, abs_err); 
    end = std::chrono::high_resolution_clock::now();

    elapsed = end-start;
    cout << elapsed.count() << " ";
    fout << elapsed.count() << " ";

    cout << endl;
    fout << endl;

  }


  return 0;
}
