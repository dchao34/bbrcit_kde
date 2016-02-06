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

void transform(double &x, double &y, 
               double mx, double my, 
               double sx, double sy, 
               double degrees) {
  x *= sx; y *= sy;
  x += mx; y += my;
  double costh = cos(degrees * M_PI / 180.);
  double sinth = sin(degrees * M_PI / 180.);
  x = costh*x - sinth*y;
  y = sinth*x + costh*y;
}

template<typename PointT>
bool ReverseExactLexicoLess(const PointT &lhs, const PointT &rhs) {
  int i = 0; while (i < lhs.dim() && lhs[lhs.dim()-i-1] == rhs[lhs.dim()-i-1]) { ++i; }
  return i != lhs.dim() && lhs[lhs.dim()-i-1] < rhs[lhs.dim()-i-1];
}

int main() {

  using KernelDensityType = KernelDensity<2>;
  using DataPointType = typename KernelDensityType::DataPointType;

  ofstream fout1("test_kde5_data.csv");
  ofstream fout2("test_kde5_kde.csv");

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
    double x = gaussian(e), y = gaussian(e), u = uniform(e);
    if (u < 0.5) {
      transform(x, y, 1, 1, 0.5, 0.3, 30);
    } else {
      transform(x, y, -1, -1, 0.5, 0.3, -30);
    }
    fout1 << x << " " << y << endl;
    data.push_back({{x,y}});
  }

  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;

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

  vector<DataPointType> queries;

  double start_x = -2, end_x = 2;
  double start_y = -2, end_y = 2;
  int x_steps = 100, y_steps = 100;
  double delta_x = (end_x - start_x) / x_steps;
  double delta_y = (end_y - start_y) / y_steps;

  for (int j = 0; j < y_steps; ++j) {

    double y_coord = start_y + j * delta_y;
    double x_coord = start_x;

    for (int i = 0; i < x_steps; ++i) {
      x_coord = start_x + i * delta_x;
      queries.push_back({{x_coord, y_coord}});
    }
  }

  for (int i = 0; i < x_steps; ++i) { fout2 << start_x + i * delta_x << " "; } fout2 << endl;
  for (int i = 0; i < y_steps; ++i) { fout2 << start_y + i * delta_y << " "; } fout2 << endl;

  start = std::chrono::high_resolution_clock::now();
  kde.eval(queries, rel_err, abs_err);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end-start;
  cout << "runtime: " << elapsed.count() << " seconds" << endl;
  cout << endl;

  std::sort(queries.begin(), queries.end(), ReverseExactLexicoLess<DataPointType>); 
  for (size_t i = 0; i < queries.size(); ++i) {
    if (i % x_steps == 0 && i) { fout2 << std::endl; }
    fout2 << queries[i].attributes().middle() << " ";
  }
  fout2 << std::endl;

  return 0;
}
