#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <vector>
#include <limits>
#include <complex>
#include <algorithm>
#include <cassert>

#include <ProdKde2d.h>
#include <fft.h>

using namespace std;

double ProdKde2d::operator()(double x1, double x2) {

  double result = 0.0;
  for (auto &p : sample) {
    result += gauss_kernel_1d((x1 - p.first)/h1) * gauss_kernel_1d((x2-p.second)/h2);
  }
  return result /= (sample.size() * h1 * h2);
}

double ProdKde2d::evaluate_marginal(double x, bool dim1) {

  decltype(get_first) *f = dim1 ? get_first : get_second;
  const double &h = dim1 ? h1 : h2;

  double result = 0.0;
  for (auto it = sample.begin(); it != sample.end(); ++it) {
    result += gauss_kernel_1d((x - f(it))/h);
  }
  return result /= (sample.size() * h);
}

double ProdKde2d::evaluate_marginal_se(double x, bool dim1) {

  decltype(get_first) *f = dim1 ? get_first : get_second;
  const double &h = dim1 ? h1 : h2;

  double mean = 0.0;
  for (auto it = sample.begin(); it != sample.end(); ++it) {
    mean += gauss_kernel_1d((x - f(it))/h);
  }
  mean /= (sample.size() * h);

  double se = 0.0;
  for (auto it = sample.begin(); it != sample.end(); ++it) {
    se += pow(gauss_kernel_1d((x - f(it))/h) - mean, 2.0);
  }
  se /= (sample.size() - 1);
  se = sqrt(se / sample.size());

  return se;
}

void ProdKde2d::compute_sample_stats() {
  mhat1 = 0; mhat2 = 0;
  for (auto &p : sample) { mhat1 += p.first; mhat2 += p.second; }
  mhat1 /= sample.size(); mhat2 /= sample.size();

  shat1 = 0; shat2 = 0;
  for (auto &p : sample) {
    double d1 = (p.first - mhat1);
    double d2 = (p.second - mhat2);
    shat1 += d1 * d1;
    shat2 += d2 * d2;
  }
  shat1 /= (sample.size() - 1);
  shat2 /= (sample.size() - 1);
}

void ProdKde2d::cv(vector<double> &results, double h, bool cv_x1) {

  decltype(get_first) *f = cv_x1 ? get_first : get_second;

  double sum = 0.0, cv = 0.0;
  auto n = sample.size();
  for (auto it1 = sample.begin(); it1 != sample.end(); ++it1) {
    for (auto it2 = sample.begin(); it2 != sample.end(); ++it2) {
      if (it1 == it2) continue;
      sum += gauss_kernel_star((f(it1) - f(it2)) / h);
    }
  }
  cv = 1/(h*n*n)*sum + 2/(n*h*sqrt(2*M_PI));

  results.clear();
  results.push_back(h);
  results.push_back(sum);
  results.push_back(cv);
}

void ProdKde2d::deduce_fft_grid_constants(
  double &a, double &b, double &delta, 
  vc_size_t M, double margin, bool cv_x1) {

  decltype(get_first) *f = cv_x1 ? get_first : get_second;

  // find lower and upper bounds on the data sample
  a = f(sample.begin()); b = a;
  for (auto it = sample.begin(); it != sample.end(); ++it) {
    double x = f(it);
    if (x < a) a = x;
    if (x > b) b = x;
  }

  a -= margin; b += margin;
  delta = (b - a) / M;
}

void ProdKde2d::discretize_data(
  vector<double> &x_disc,
  double a, double b, double delta, 
  vc_size_t M, vc_size_t N, bool cv_x1) {

  decltype(get_first) *f = cv_x1 ? get_first : get_second;

  // create bins
  vector<double> t(M);
  for (vc_size_t k = 0; k < M; ++k) {
    t[k] = a + k * delta;
  }

  // fill the histogram 
  x_disc.clear(); x_disc = vector<double>(M, 0.0);
  for (auto is = sample.begin(); is != sample.end(); ++is) {
    double x = f(is);
    auto it = lower_bound(t.begin(), t.end(), x);
    if (it != t.end()) {
      assert(it != t.begin());
      vc_size_t k = static_cast<vc_size_t>(it - t.begin()) - 1;
      x_disc[k] += (t[k+1] - x)/(N*delta*delta);
      x_disc[k+1] += (x - t[k])/(N*delta*delta);
    } else {
      x_disc[M-1] += (b - x)/(N*delta*delta);
    }
  }
}

// In Silverman's notation:
// Input: discretized data x_disc[k], k=0,...,M-1.
// Ouput: Y_l's, where Y_l = y[M/2+l], l=-M/2,...,M/2-1.
void ProdKde2d::fft_forward(
    const vector<double> &x_disc,
    vector<cplex> &y, unsigned r) {

  vc_size_t M = x_disc.size();

  // preprocess x_disc by shifting. this is necessary since the fft
  // routine used is in a different convention (CLRS). 
  vector<cplex> x_cplex(M);
  for (vector<double>::size_type k = 0; k < M; ++k) {
    double sign = 1.0;
    if (k % 2) sign *= -1.0;
    x_cplex[k] = cplex(sign*x_disc[k], 0.0);
  }

  y.clear(); y = vector<cplex>(M);
  fft(x_cplex, y, r);
  for (auto &c : y) { c /= M; }
}

void ProdKde2d::grid_evaluate_marginal(
    vector<pair<double,double>> &results, 
    bool dim1, unsigned r) {

  const double &h = dim1 ? h1 : h2;

  // cache constants 
  vc_size_t M = 0x1 << r; 
  vp_size_t N = sample.size();

  // choice of a and b are important. may need to let the user set `margin`
  double a, b, delta;
  deduce_fft_grid_constants(a, b, delta, M, 3*h, dim1);

  // discretize data 
  vector<double> x_disc;
  discretize_data(x_disc, a, b, delta, M, N, dim1);

  // fft
  vector<cplex> y;
  fft_forward(x_disc, y, r);

  // ifft
  for (vc_size_t l = 1; l <= M/2; ++l) {
    double s = 2*M_PI*l/(b-a);
    double t = exp(-0.5*h*h*s*s);
    y[M/2-l] *= t;
    if (l < M/2) 
      y[M/2+l] *= t;
  }
  vector<cplex> f(M);
  ifft(y, f, r);

  // save results
  results.clear(); results = vector<pair<double, double>>(M);
  for (vc_size_t k = 0; k < M; ++k) {
    double sign = 1.0;
    if (k % 2) sign *= -1.0;
    f[k] *= sign * M;
    results[k] = { a + k*delta, f[k].real() };
  }

}

// TODO: score of fcv and cv differ by ~1
void ProdKde2d::fcv(vector<double> &results, double h, unsigned r, bool cv_x1) {

  // cache constants 
  vc_size_t M = 0x1 << r; 
  vp_size_t N = sample.size();

  // choice of a and b are important. may need to let the user set `margin`
  double a, b, delta;
  deduce_fft_grid_constants(a, b, delta, M, 10.0, cv_x1);

  // discretize data 
  vector<double> x_disc;
  discretize_data(x_disc, a, b, delta, M, N, cv_x1);

  // fft
  vector<cplex> y;
  fft_forward(x_disc, y, r);

  // compute cv score
  double sum = 0.0;
  for (vc_size_t l = 1; l <= M/2; ++l) {
    double s = 2*M_PI*l/(b-a);
    double t = exp(-0.5*h*h*s*s);
    sum += (t*t-2*t)*norm(y[M/2-l]);
  }
  double cv = 2* (b-a) * sum - 1;
  cv += 2/(N*h*sqrt(2*M_PI));

  results.clear();
  results.push_back(h);
  results.push_back(sum);
  results.push_back(cv);

}
