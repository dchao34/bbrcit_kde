#include <iostream>
#include <fft.h>

using namespace std;

using cplex = complex<double>;
using vc_size_t = vector<cplex>::size_type;

void print(const vector<cplex> &a, bool real=false) {
  cout << "(";
  for (vc_size_t i = 0; i < a.size(); ++i) {
    if (i != 0) cout << ", ";
    if (real) {
      cout << a[i].real();
    } else {
      cout << a[i];
    }
  }
  cout << ")";
}

int main() {

  cout << "Multiplying polynomial (1.0, 3.0, 3.0, 1.0) by itself: " << endl;
  vector<cplex> a = { 1.0, 3.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
  vector<cplex> b = a;
  vector<cplex> A(8);
  vector<cplex> B(8);
  vector<cplex> C(8);
  fft(a, A, 3);
  fft(a, B, 3);
  for (unsigned i = 0; i < 8; i++) {
    C[i] = A[i] * B[i];
  }

  vector<cplex> c(8);
  ifft(C, c, 3);
  print(c, true); cout << endl;
  cout << endl;

  cout << "a = (1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0): " << endl;
  a = { 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0 };
  vector<cplex> y(a.size());
  fft(a, y, 3);
  cout << "DFT(a) = ";
  print(y); cout << endl;
  cout << "IDFT(DFT(a)) = ";
  ifft(y, a, 3);
  print(a); cout << endl;
  cout << endl;

  cout << "a = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0): " << endl;
  a = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
  y = a;
  fft(a, y, 3);
  cout << "DFT(a) = ";
  print(y); cout << endl;
  cout << "IDFT(DFT(a)) = ";
  ifft(y, a, 3);
  print(a); cout << endl;
  cout << endl;

  cout << "a = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0): " << endl;
  a = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
  y = a;
  fft(a, y, 4);
  cout << "DFT(a) = ";
  print(y); cout << endl;
  cout << "IDFT(DFT(a)) = ";
  ifft(y, a, 4);
  print(a); cout << endl;
  cout << endl;

  return 0;
}
