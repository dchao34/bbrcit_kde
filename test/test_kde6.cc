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

  using EpanKdeType = KernelDensity<2,double,EpanechnikovKernel<2,double>>;
  using GaussKdeType = KernelDensity<2,double,GaussianKernel<2,double>>;
  using DataPointType = typename EpanKdeType::DataPointType;

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed;

  mt19937 e;
  normal_distribution<> gaussian(0, 1);

  ofstream fout("test_kde6.csv");

  vector<DataPointType> data, queries;
  double rel_err = 1e-6; double abs_err = 1e-10;

  int n_samples;
  for (int k = 7; k < 22; ++k) {

    cout << k << " ";
    fout << k << " ";

    data.clear();

    n_samples = std::pow(2, k);
    for (int i = 0; i < n_samples; ++i) { 
      data.push_back({{gaussian(e), gaussian(e)}}); 
    }

    // Epanechnikov 
    // ------------

    // kde build
    start = std::chrono::high_resolution_clock::now();
    EpanKdeType *epan_kde = new EpanKdeType(data); 
    end = std::chrono::high_resolution_clock::now();
    elapsed = end-start;
    cout << elapsed.count() << " ";
    fout << elapsed.count() << " ";


    // compute optimal bandwidth
    double epan_bw = 0.01;
    //epan_bw = 2.345 * std::pow(n_samples, -0.2);
    cout << epan_bw << " ";
    fout << epan_bw << " ";
    epan_kde->kernel().set_bandwidth(epan_bw);

    // direct evaluation
    if (k < 13) {
      queries = data;
      start = std::chrono::high_resolution_clock::now();
      for (auto &p : queries) { 
        auto result = epan_kde->direct_eval(p); 
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

    // dual tree evaluation
    queries = data;

    start = std::chrono::high_resolution_clock::now();
    epan_kde->eval(queries, rel_err, abs_err); 
    end = std::chrono::high_resolution_clock::now();

    elapsed = end-start;
    cout << elapsed.count() << " ";
    fout << elapsed.count() << " ";

    delete epan_kde;

    // Gaussian 
    // --------

    // kde build
    GaussKdeType *gauss_kde = new GaussKdeType(data); 

    // compute optimal bandwidth
    double gauss_bw = epan_bw / 2.345;
    //gauss_bw = 1.059 * std::pow(n_samples, -0.2);
    cout << gauss_bw << " ";
    fout << gauss_bw << " ";
    gauss_kde->kernel().set_bandwidth(epan_bw);

    // direct evaluation
    if (k < 13) {
      queries = data;
      start = std::chrono::high_resolution_clock::now();
      for (auto &p : queries) { 
        auto result = gauss_kde->direct_eval(p); 
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
      auto result = gauss_kde->eval(p, rel_err, abs_err); 
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
    gauss_kde->eval(queries, rel_err, abs_err); 
    end = std::chrono::high_resolution_clock::now();

    elapsed = end-start;
    cout << elapsed.count() << " ";
    fout << elapsed.count() << " ";


    delete gauss_kde;

    cout << endl;
    fout << endl;

  }


  return 0;
}
