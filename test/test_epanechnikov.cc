#include <iostream>
#include <fstream>
#include <vector>

#include <EpanechnikovKernel.h>
#include <Point.h>

using namespace std;
using bbrcit::EpanechnikovKernel;
using Point1d = bbrcit::Point<1>;

vector<Point1d> generate_grid1d(double start, double end, int steps) {
  vector<Point1d> result;
  double delta = (end - start) / steps;
  for (int i = 0; i < steps; ++i) {
    result.push_back(Point1d({start+i*delta}));
  }
  return result;
}

int main() {
  EpanechnikovKernel<1> k;
  cout << k.dim() << " " << k.normalization << endl;

  ofstream fout("test_epanechnikov1d.csv");
  vector<Point1d> grid = generate_grid1d(-2, 2, 1000);
  for (const auto &p : grid) {
    fout << p[0] << " " << k.eval(p) << endl;
  }
  
  return 0;
}
