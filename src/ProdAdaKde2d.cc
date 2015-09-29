#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include <cmath>

#include <ProdAdaKde2d.h>

using namespace std;

ProdAdaKde2d::ProdAdaKde2d(
    string fname, 
    double bw1, double bw2, 
    double a1, double a2, unsigned r) : 
  ProdKde2d(fname, bw1, bw2), 
  alpha1(a1), alpha2(a2) {

  vector<double> bwidths1;
  compute_adaptive_bandwidths(bwidths1, true, r);
  vector<double> bwidths2;
  compute_adaptive_bandwidths(bwidths2, false, r);
  assert(bwidths1.size() == bwidths2.size());

  bwidths.clear(); bwidths.reserve(bwidths1.size());
  for (vp_size_t i = 0; i < bwidths1.size(); ++i) {
    bwidths.push_back({ bwidths1[i], bwidths2[i] });
  }

}

void ProdAdaKde2d::compute_adaptive_bandwidths(
    vector<double> &results, bool dim1, unsigned r) {

  decltype(get_first) *f = dim1 ? get_first : get_second;
  double alpha = dim1 ? alpha1 : alpha2;
  
  vector<pair<double, double>> grid_estimate;
  grid_evaluate_marginal(grid_estimate, dim1, r);

  vector<double> pilot_estimate;
  for (auto is = sample.begin(); is != sample.end(); ++is) {

    pair<double, double> x = { f(is), 0.0 };
    auto it = lower_bound(grid_estimate.begin(), grid_estimate.end(), x);

    double estimate; auto prev = it - 1;
    if (it != grid_estimate.end()) {
      assert(it != grid_estimate.begin());
      estimate = prev->second + 
                 (it->second - prev->second) * (f(is) - prev->first) / (it->first - prev->first);
    } else {
      estimate = prev->second;
    }
    assert(estimate != 0.0);

    pilot_estimate.push_back(estimate);
  }

  double g = 0.0;
  for (auto e : pilot_estimate) {
    g += log(e);
  }
  g /= sample.size();
  g = exp(g);

  results.clear(); results.reserve(sample.size());
  for (auto e : pilot_estimate) {
    results.push_back(pow(e/g, -alpha));
  }

}

// evaluate the density estimate at a new point
double ProdAdaKde2d::operator()(double x1, double x2) {
  double result = 0.0;
  for (vp_size_t i = 0; i < sample.size(); ++i) {
    double h1a = h1 * bwidths[i].first; double h2a = h2 * bwidths[i].second;
    result += 1/(h1a*h2a) * gauss_kernel_1d((x1-sample[i].first)  / h1a) *
                            gauss_kernel_1d((x2-sample[i].second) / h2a);
  }
  return result /= sample.size();
}

double ProdAdaKde2d::evaluate_marginal(double x, bool dim1) {

  decltype(get_first_a) *f = dim1 ? get_first_a : get_second_a;
  double h = dim1 ? h1 : h2;

  double result = 0.0;
  for (vp_size_t i = 0; i < sample.size(); ++i) {
    double ha = h * f(bwidths[i]);
    result += 1/ha * gauss_kernel_1d((x-f(sample[i])) / ha);
  }
  return result /= sample.size();
}

void ProdAdaKde2d::cv(vector<double> &results, double h, bool dim1) {

  decltype(get_first_a) *f = dim1 ? get_first_a : get_second_a;

  double sum1 = 0.0;
  for (vp_size_t i = 0; i < sample.size(); ++i) {
    for (vp_size_t j = 0; j < sample.size(); ++j) {
      if (i == j) continue;
      double ha = h * f(bwidths[j]);
      sum1 += gauss_kernel_1d((f(sample[i]) - f(sample[j])) / ha) / ha;
    }
  }
  sum1 /= sample.size() * (sample.size() - 1);

  double sum2 = 0.0;
  for (vp_size_t i = 0; i < sample.size(); ++i) {
    for (vp_size_t j = 0; j < sample.size(); ++j) {
      double ha = h * f(bwidths[i]);
      double hb = h * f(bwidths[j]);
      double s = sqrt(ha*ha + hb*hb);
      sum2 += gauss_kernel_1d(f(sample[i]) - f(sample[j]), s);
    }
  }
  sum2 /= (sample.size() * sample.size());

  double cv = sum2 - 2 * sum1;

  results.clear();
  results.push_back(h);
  results.push_back(sum1);
  results.push_back(sum2);
  results.push_back(cv);
}
