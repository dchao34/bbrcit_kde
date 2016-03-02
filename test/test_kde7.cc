#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>

#include <KernelDensity.h>
#include <Kernels/EpanechnikovKernel.h>
#include <Point.h>

#include "kde_test_utils.h"

using namespace std;

using bbrcit::KernelDensity;
using bbrcit::EpanechnikovKernel;
using bbrcit::Point;

int main() {

  using EpanKdeType = KernelDensity<2,EpanechnikovKernel<2>>;
  using DataPointType = typename EpanKdeType::DataPointType;

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed;

  mt19937 e;
  normal_distribution<> gaussian(0, 1);

  ofstream fout("test_kde7.csv");

  vector<DataPointType> data, queries;
  double rel_err = 1e-6; double abs_err = 1e-10;

  vector<double> bandwidths; bandwidths.push_back(0.1); 
  for (int i = 1; i < 6; ++i) {
    bandwidths.push_back(bandwidths[i-1]*0.1);
  }

  int n_samples;
  for (int k = 7; k < 16; ++k) {

    cout << k << " ";
    fout << k << " ";

    n_samples = std::pow(2, k);
    for (int i = 0; i < n_samples; ++i) { 
      data.push_back({{gaussian(e), gaussian(e)}}); 
    }

    // kde build
    EpanKdeType *epan_kde = new EpanKdeType(data); 

    for (auto bw : bandwidths) {

      epan_kde->kernel().set_bandwidth(bw);

      // naive evaluation
      /*
      queries = data;
      start = std::chrono::high_resolution_clock::now();
      for (auto &p : queries) { 
        auto result = epan_kde->naive_eval(p); 
        auto attr = p.attributes();
        attr.set_lower(result);
        attr.set_upper(result);
        p.set_attributes(attr);
      }
      end = std::chrono::high_resolution_clock::now();

      elapsed = end-start;
      cout << elapsed.count() << " ";
      fout << elapsed.count() << " ";

      // single tree evaluation
      queries = data;
      start = std::chrono::high_resolution_clock::now();
      for (auto &p : queries) { 
        auto result = epan_kde->eval(p, rel_err, abs_err); 
        auto attr = p.attributes();
        attr.set_lower(result);
        attr.set_upper(result);
        p.set_attributes(attr);
      }
      end = std::chrono::high_resolution_clock::now();

      elapsed = end-start;
      cout << elapsed.count() << " ";
      fout << elapsed.count() << " ";
      */


      // dual tree evaluation
      queries = data;

      start = std::chrono::high_resolution_clock::now();
      epan_kde->eval(queries, rel_err, abs_err); 
      end = std::chrono::high_resolution_clock::now();

      elapsed = end-start;
      cout << elapsed.count() << " ";
      fout << elapsed.count() << " ";


    }

    delete epan_kde;

    cout << endl;
    fout << endl;

  }


  return 0;
}
