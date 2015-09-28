#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <vector>
#include <complex>
#include <cassert>

using namespace std;

using cplex = complex<double>;
using vc_size_t = vector<cplex>::size_type;

static cplex I = cplex(0.0, 1.0);

vc_size_t bit_reverse(vc_size_t n, unsigned r) {
  for (unsigned i = 0; i < r/2; ++i) {
    vc_size_t lo = n >> i & 0x1ULL;
    vc_size_t hi = n >> (r-1-i) & 0x1ULL;
    if (hi^lo) {
      n ^= 0x1ULL << i | 0x1ULL << (r-1-i);
    }
  }
  return n;
}

void bit_reverse_copy(const vector<cplex> &a, vector<cplex> &y, unsigned r) {
  assert(a.size() == y.size());
  for (vc_size_t i = 0; i < a.size(); ++i) {
    y[i] = a[bit_reverse(i, r)];
  }
}

void fft(const vector<cplex> &a, vector<cplex> &y, unsigned r) {
  vc_size_t n = a.size();
  assert(n == (0x1ull << r));
  assert(a.size() == y.size());

  bit_reverse_copy(a, y, r);
  for (unsigned s = 1; s <= r; ++s) {
    vc_size_t m = 0x1ull << s;
    cplex w_p = exp(2*M_PI/m*I);
    for (vc_size_t k = 0; k < n; k += m) {
      cplex w = cplex(1.0, 0.0);
      for (vc_size_t j = 0; j < m/2; ++j) {
        cplex t = w * y[k+j+m/2];
        cplex u = y[k+j];
        y[k+j] = u + t;
        y[k+j+m/2] = u - t;
        w *= w_p;
      }
    }
  }
}

void ifft(const vector<cplex> &a, vector<cplex> &y, unsigned r) {
  vc_size_t n = a.size();
  assert(n == (0x1ull << r));
  assert(a.size() == y.size());

  bit_reverse_copy(a, y, r);
  for (unsigned s = 1; s <= r; ++s) {
    vc_size_t m = 0x1 << s;
    cplex w_p = exp(-2*M_PI/m*I);
    for (vc_size_t k = 0; k < n; k += m) {
      cplex w = cplex(1.0, 0.0);
      for (vc_size_t j = 0; j < m/2; ++j) {
        cplex t = w * y[k+j+m/2];
        cplex u = y[k+j];
        y[k+j] = u + t;
        y[k+j+m/2] = u - t;
        w *= w_p;
      }
    }
  }
  for (auto &c : y) {
    c /= n;
  }
}
