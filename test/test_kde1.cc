#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

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

int main() {

  using KernelDensityType = KernelDensity<2>;
  using DataPointType = typename KernelDensityType::DataPointType;

  ofstream fout1("test_kde1_data.csv");
  ofstream fout2("test_kde1_kde.csv");

  cout << "generating data" << endl;
  default_random_engine e;
  normal_distribution<> gaussian(0.0, 1.0);
  uniform_real_distribution<> uniform(0.0, 1.0);
  vector<DataPointType> data;
  for (int i = 0; i < 10000; ++i) {
    double x = gaussian(e), y = gaussian(e), u = uniform(e);
    if (u < 0.5) {
      transform(x, y, 1, 1, 0.5, 0.3, 30);
    } else {
      transform(x, y, -1, -1, 0.5, 0.3, -30);
    }
    fout1 << x << " " << y << endl;
    data.push_back({{x,y}});
  }

  KernelDensityType kde(data, 0.5, 64);
  //kde.eval({0,0}, 0);

  cout << "generating grid/evaluating kde" << endl;
  double start_x = -2, end_x = 2;
  double start_y = -2, end_y = 2;
  int x_steps = 100, y_steps = 100;
  double delta_x = (end_x - start_x) / x_steps;
  double delta_y = (end_y - start_y) / y_steps;

  for (int i = 0; i < x_steps; ++i) { fout2 << start_x + i * delta_x << " "; }
  fout2 << endl;
  for (int i = 0; i < y_steps; ++i) { fout2 << start_y + i * delta_y << " "; }
  fout2 << endl;

  for (int j = 0; j < y_steps; ++j) {

    double y_coord = start_y + j * delta_y;
    double x_coord = start_x;

    fout2 << kde.naive_eval({x_coord, y_coord}); 
    for (int i = 1; i < x_steps; ++i) {
      x_coord = start_x + i * delta_x;
      //fout2 << " " << kde.naive_eval({x_coord, y_coord}); 
      fout2 << " " << kde.eval({x_coord, y_coord}, 1e-12); 
    }
    fout2 << endl;
  }

  return 0;
}
