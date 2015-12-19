#include <iostream>
#include <fstream>

#include <Kernel.h>
#include <Point.h>

using namespace std;

int main() {

  GaussKernel1d k;

  cout << "Test 1: Single point evaluations. " << endl;
  Point1d p; cout << k(p) << " (should be close to 0.398942) " << endl;
  p = 1.5; cout << k(p) << " (should be close to 0.129518) " << endl;

  cout << endl;

  string fname("gauss_kernel_trace.csv");
  ofstream fout(fname);
  int n = 1000; double bound = 5.0;
  for (int i = 0; i < n; ++i) {
    p = -bound + i * 2 * bound / n;
    fout << p[0] << " " << k(p) << endl;
  }
  cout << "Test 2: Kernel trace saved to " << fname << endl;

  return 0;
}
