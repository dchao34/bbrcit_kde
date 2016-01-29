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
  int n_data = 1000;
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
  //KernelDensityType kde(data, 0.01, 2);

  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;
  

  // evaluate kde at grid points
  cout << "evaluating kde" << endl;
  start = std::chrono::high_resolution_clock::now();

  vector<size_t> leaves_visited; leaves_visited.reserve(data.size());
  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < data.size(); ++i) {
    size_t cnt = 0;
    kde.eval(data[i].point(), 1e-3, 1e-1, cnt); 
    leaves_visited.push_back(cnt);
  }

  end = std::chrono::high_resolution_clock::now();
  cout << "min/max leaves visited: "; 
  auto it = std::min_element(leaves_visited.begin(), leaves_visited.end());
  cout << (*it) << ", ";
  it = std::max_element(leaves_visited.begin(), leaves_visited.end());
  cout << (*it) << endl;
  cout << "mean leaves visited: "; 
  cout << accumulate(leaves_visited.begin(), leaves_visited.end(), 0.0) / leaves_visited.size() << endl;
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;


  return 0;
}
