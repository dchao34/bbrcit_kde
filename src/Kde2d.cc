#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <iostream>
#include <numeric>
#include <random>
#include <algorithm>
#include <cassert>

#include "gauss_legendre.h"
#include "file_io_utils.h"
#include "Kde2d.h"

using namespace std;

Kde2d::Kde2d(string data_fname, double bw1, double bw2) {

  h1 = bw1; h2 = bw2;

  ifstream fin;
  open_for_reading(fin, data_fname);

  // read the file line by line
  string line; double x1, x2;
  while (getline(fin, line)) {

    // read each line column by column
    istringstream sin(line);
    sin >> x1 >> x2;

    sample.push_back({x1, x2});
  }

  fin.close();
}

double Kde2d::operator()(double x1, double x2) {
  double result = 0.0;
  for (auto &p : sample) {
    result += gauss2d(x1, x2, p.first, p.second);
  }
  return result / sample.size();
}

// evaluate the density at (x1, x2) using all sample pts other 
// than those in [bid, eidx).
double Kde2d::excluded_eval(double x1, double x2, sample_no bidx, sample_no eidx) {

  assert(bidx <= eidx);

  double result = 0.0;
  for (sample_no i = 0; i < bidx; i++) {
    result += gauss2d(x1, x2, sample[i].first, sample[i].second);
  }
  for (sample_no i = eidx; i < sample.size(); i++) {
    result += gauss2d(x1, x2, sample[i].first, sample[i].second);
  }
  return result / (sample.size() - (eidx - bidx));
}

double f2(double x1, double x2, void *kde_obj_addr) {
  Kde2d *kde = (Kde2d *) kde_obj_addr;
  return kde->eval2(x1, x2);
}

void Kde2d::cv(
    ostream &os, 
    const vector<pair<double, double>> &candidates, 
    double x_low, double x_high, double y_low, double y_high,
    int qgauss_n) {

  // permute the data
  random_device rd;
  mt19937_64 e(rd());
  shuffle(sample.begin(), sample.end(), e);

  for (auto &h : candidates) {

    set_h1(h.first); set_h2(h.second);

    double s1 = gauss_legendre_2D_cube(qgauss_n, f2, this, x_low, x_high, y_low, y_high);

    double s2 = 0;
    sample_no n = sample.size();
    for (sample_no i = 0; i < n; i++) {

      auto bidx = i; auto eidx = i+1;

      for (sample_no i = bidx; i < eidx; i++) {
        s2 += excluded_eval(sample[i].first, sample[i].second, bidx, eidx);
      }
    }
    s2 *= (2.0 / n);

    os << "h1 = " << h1 << ", h2 = " << h2;
    os << " (" << s1 << ", " << s2 << ", " << s1 - s2 << ")" << endl;
  }

}
