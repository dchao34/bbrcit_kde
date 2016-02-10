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

// this test is to check that the single and dual tree algorithms
// get the same answers, at least visually, as compared to the 
// naive algorithm. 
//
// input: 2d data points sampled from a bimodal gaussian, all 
//        equally weighted. 
// kernel: 2d epanechnikov.

int main() {

  using KernelDensityType = KernelDensity<2>;
  using DataPointType = typename KernelDensityType::DataPointType;

  ofstream fout;

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed;

  // 1. generate the data
  int n_samples = 10000;
  cout << "+ generating " << n_samples << " points. " << endl;

  default_random_engine e;
  vector<DataPointType> data;

  start = std::chrono::high_resolution_clock::now();
  generate_bimodal_gaussian(e, data, n_samples, 
                            1, 1, 0.5, 0.3, 30, 
                            -1, -1, 0.5, 0.3, -30);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "  runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;

  fout.open("test_kde3_data.csv");
  write_scatter_data(fout, data);
  fout.close();

  // 2. build the kernel density estimator
  cout << "+ building the kde (based on a kdtree)" << endl;

  start = std::chrono::high_resolution_clock::now();
  KernelDensityType kde(data); 
  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "  runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;
  
  // configure the kernel
  kde.kernel().set_bandwidth(0.2);
  

  // 3. evaluate the kde fixed grid points
  double rel_err = 1e-6; double abs_err = 1e-10;
  cout << "+ evaluating the kde with rel_err = " << rel_err; 
  cout << ", abs_err = " << abs_err << endl;
  cout << endl;

  // generate the grid
  vector<DataPointType> grid, queries;
  double start_x = -2, end_x = 2; int steps_x = 100;
  double start_y = -2, end_y = 2; int steps_y = 100;
  generate_2dgrid(grid, start_x, end_x, steps_x, start_y, end_y, steps_y);

  // naive evaluation
  cout << "  naive"; 
  queries = grid;

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
  cout << " runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;

  fout.open("test_kde3_naive.csv");
  write_kde2d_result(fout, queries, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);
  fout.close();

  // single tree evaluation
  cout << "  single tree"; 
  queries = grid;

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
  cout << " runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;

  fout.open("test_kde3_single.csv");
  write_kde2d_result(fout, queries, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);
  fout.close();

  // dual tree evaluation
  cout << "  dual tree"; 
  queries = grid;

  start = std::chrono::high_resolution_clock::now();
  kde.eval(queries, rel_err, abs_err); 
  end = std::chrono::high_resolution_clock::now();

  elapsed = end-start;
  cout << " runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;

  fout.open("test_kde3_dual.csv");
  write_kde2d_result(fout, queries, 
                     start_x, end_x, steps_x,
                     start_y, end_y, steps_y);
  fout.close();

  return 0;
}
