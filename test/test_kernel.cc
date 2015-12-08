#include <Kernel.h>
#include <Point.h>
#include <iostream>

using namespace std;

int main() {

  GaussKernel1d k;

  // 0.398942
  Point1d p;
  cout << k(p) << endl;

  // 0.129518
  p[0] = 1.5;
  cout << k(p) << endl;

  return 0;
}
