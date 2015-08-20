#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <vector>
#include <limits>
#include <complex>
#include <algorithm>
#include <cassert>

#include "ProdKde2d.h"
#include "fft.h"

using namespace std;

using vp_it = vector<pair<double,double>>::iterator;
double get_first(vp_it p) { return p->first; }
double get_second(vp_it p) { return p->second; }

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
  cv = 1/(h*n*n)*sum + 2/(n*h)/sqrt(2*M_PI);

  results.clear();
  results.push_back(h);
  results.push_back(sum);
  results.push_back(cv);
}

void ProdKde2d::cv(ostream &os, const vector<double> &candidates, bool cv_h1) {

  decltype(get_first) *f = cv_h1 ? get_first : get_second;
  if (cv_h1) {
    os << "cross validating h1:" << endl;
  } else {
    os << "cross validating h2:" << endl;
  }

  double cv_min = 0.0, cv_argmin = 0.0;
  for (auto &h : candidates) {
    double sum = 0.0, cv = 0.0;
    auto n = sample.size();
    for (auto it1 = sample.begin(); it1 != sample.end(); ++it1) {
      for (auto it2 = sample.begin(); it2 != sample.end(); ++it2) {
        sum += gauss_kernel_star((f(it1) - f(it2)) / h);
      }
    }

    cv = 1/(h*n*n)*sum + 2/(n*h)/sqrt(2*M_PI);
    if (cv < cv_min) { cv_min = cv; cv_argmin = h; }

    os << "h = " << h;
    os << ", sum = " << sum;
    os << ", cv = " << cv << endl;
  }
  os << "best cv = " << cv_min << ", best h = " << cv_argmin << endl;
}

void ProdKde2d::fcv(vector<double> &results, double h, unsigned r, bool cv_x1) {

  using cplex = complex<double>;
  using vc_size_t = vector<cplex>::size_type;
  using vp_size_t = vector<pair<double,double>>::size_type;

  decltype(get_first) *f = cv_x1 ? get_first : get_second;

  // cache constants 
  vc_size_t M = 0x1 << r;
  vp_size_t N = sample.size();

  double a, b;
  a = f(sample.begin()); b = a;
  for (auto it = sample.begin(); it != sample.end(); ++it) {
    double x = f(it);
    if (x < a) a = x;
    if (x > b) b = x;
  }
  a -= 4*h; b += 4*h;

  double delta = (b - a) / M;

  // discretize data 
  vector<double> t(M);
  for (vc_size_t k = 0; k < M; ++k) {
    t[k] = a + k * delta;
  }

  vector<double> x_disc(M);
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
  for (vector<double>::size_type k = 0; k < x_disc.size(); ++k) {
    if (k % 2) x_disc[k] *= -1;
  }

  vector<cplex> x_cplex(M);
  for (vc_size_t i = 0; i < M; ++i) {
    x_cplex[i] = cplex(x_disc[i], 0.0);
  }

  // fft
  vector<cplex> y(M);
  fft(x_cplex, y, r);
  for (auto &c : y) { c /= M; }

  // compute cv score
  double sum = 0.0;
  for (vc_size_t l = 1; l <= M/2; ++l) {
    double s = 2*M_PI*l/(b-a);
    double t = exp(-0.5*h*h*s*s);
    sum += (t*t-2*t)*norm(y[M/2-l]);
  }
  double cv = (b-a) * sum + 1/(N*h*sqrt(2*M_PI));
  cv = 2 * cv -1;

  results.clear();
  results.push_back(h);
  results.push_back(sum);
  results.push_back(cv);

}
