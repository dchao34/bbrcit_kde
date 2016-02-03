#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>

#include <KernelDensity.h>
#include <Point.h>

using namespace std;

using bbrcit::KernelDensity;
using bbrcit::Point;

int main() {

  using KernelDensityType = KernelDensity<2>;
  using DataPointType = typename KernelDensityType::DataPointType;

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed;

  // generating data
  int n_data = 10000;
  cout << "generating data: " << n_data << endl;
  start = std::chrono::high_resolution_clock::now();

  default_random_engine e;
  normal_distribution<> gaussian(0.0, 1.0);
  uniform_real_distribution<> uniform(0.0, 1.0);
  vector<DataPointType> data;
  for (int i = 0; i < n_data; ++i) {
    double x = gaussian(e), y = gaussian(e);
    data.push_back({{x,y}});
  }

  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;

  // building tree
  cout << "building kdtree" << endl;
  start = std::chrono::high_resolution_clock::now();

  KernelDensityType kde(data, 2.4*std::pow(n_data, -1/6.0), 2);

  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;
  

  // evaluate kde at grid points
  cout << "evaluating kde" << endl;

  start = std::chrono::high_resolution_clock::now();
  kde.eval(data, 1e-12);

  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;


  return 0;
}
