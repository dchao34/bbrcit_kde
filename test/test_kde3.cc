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

#include "kde_test_utils.h"

using namespace std;

using bbrcit::KernelDensity;
using bbrcit::Point;

int main() {

  using KernelDensityType = KernelDensity<2>;
  using DataPointType = typename KernelDensityType::DataPointType;

  ofstream fout1("test_kde3_data.csv");
  ofstream fout2("test_kde3_kde.csv");

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed;

  // generate data
  int n_samples = 100000;
  default_random_engine e;
  vector<DataPointType> data;
  cout << "generating data: " << n_samples << endl;

  start = std::chrono::high_resolution_clock::now();
  generate_bimodal_gaussian(e, data, n_samples, 
                            1, 1, 0.5, 0.3, 30, 
                            -1, -1, 0.5, 0.3, -30);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;

  write_kde2d_data(fout1, data);

  // building tree
  cout << "building kdtree" << endl;

  start = std::chrono::high_resolution_clock::now();
  KernelDensityType kde(data, 0.2, 2);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;
  

  // evaluate kde at grid points
  double rel_err = 1e-6; double abs_err = 1e-10;
  cout << "evaluating kde. rel_err = " << rel_err; 
  cout << ", abs_err = " << abs_err << endl;

  vector<DataPointType> grid;
  double start_x = -2, end_x = 2; int steps_x = 100;
  double start_y = -2, end_y = 2; int steps_y = 100;
  generate_2dgrid(grid, start_x, end_x, steps_x, start_y, end_y, steps_y);

  start = std::chrono::high_resolution_clock::now();
  for (auto &p : grid) { 
    auto result = kde.eval(p.point(), rel_err, abs_err); 
    auto attr = p.attributes();
    attr.set_lower(result);
    attr.set_upper(result);
    p.set_attributes(attr);
  }
  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;

  write_kde2d_result(fout2, grid, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);

  return 0;
}
