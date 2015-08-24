#include <iostream>
#include <iomanip>

#include <file_io_utils.h>

#include "ProdAdaKde2d.h"

using namespace std;

// evaluates kde for a given sample using specified bandwidths and grid points. 
// 2d output: 
//   + output will consists of m+2 lines with columns separated by ' '. 
//     suppose:
//     *f* is the kde
//     *m* is the total number of grid points over the 1st dimension.
//     *n* is the total number of grid points over the 2nd dimension.
//
//
//     line 1: b1 b1+s1 b1+2*s1 ... b1+(m-1)*s1
//     line 2: b2 b2+s2 b1+2*s2 ... b1+(n-1)*s2
//     line 3: f(b1, b2) f(b1+s1, b2) ... f(b1+(n-1)*s1, b2)
//     line 4: f(b1, b2+s2) f(b1+s1, b2+s2) ... f(b1+(n-1)*s1, b2+s2)
//     ...
//     line m+2: f(b1, b2+(m-1)*s2) f(b1+s1, b2+(m-1)*s2) ... f(b1+(n-1)*s1, b2+(m-1)*s2)
//
// 1d output:
//   + suppose:
//     *h1(h2)* is the bandwidth for the 1st(2nd) dimension. 
//     *a1(a2)* is the sensitivity parameter for the 1st(2nd) dimension. 
//     *f1* is the marginal kde for the 1st dimension.
//     *f2* is the marginal kde for the 2nd dimension.
//     *m* is the total number of grid points over the 1st dimension.
//     *n* is the total number of grid points over the 2nd dimension.
//     line 1: h1 h2 a1 a2
//     line 2: b1 b1+s1 b1+2*s1 ... b1+(m-1)*s1
//     line 3: b2 b2+s2 b1+2*s2 ... b1+(n-1)*s2
//     line 4: f1(b1) f1(b1+s1) ... f1(b1+(n-1)*s1)
//     line 5: f2(b2) f2(b2+s2) ... f2(b2+(n-2)*s2)
int main(int argc, char *argv[]) {

  if (argc < 7) {
    cerr << endl;
    cerr << "usage: ./kde_scan 2d_out 1d_out data_in bw1 bw2 a1 a2 r b1 e1 s1 b2 e2 s2" << endl;
    cerr << endl;
    cerr << left << setfill('.');
    cerr << setw(30) << "1. 2d_out" <<  "file to save 2d scan. " << endl;
    cerr << setw(30) << "2. 1d_out" << "file to save 1d scans. " << endl;
    cerr << setw(30) << "3. data_in" << "file containing the sample points. " << endl;
    cerr << setw(30) << "4-5. bw1 bw2" << "bandwidths in 1st and 2nd dimensions. " << endl;
    cerr << setw(30) << "6-7. a1 a2" << "adaptive sensitivity. " << endl;
    cerr << setw(30) << "8. r" << "bits in fft. (sets precision)" << endl;
    cerr << setw(30) << "9-11. b1, e1, s1" << "grid points in 1st dimension. " << endl;
    cerr << setw(30) << setfill(' ') << "" << "covers the points [<b1>, <e1>) "
                                              "in *increments* of <s1>." << endl; 
    cerr << setfill('.');
    cerr << setw(30) << "12-14. b2, e2, s2" << "grid points in 2nd dimension. " << endl;
    cerr << endl;
    return 1;
  }

  // load kde
  double h1 = atof(argv[4]), h2 = atof(argv[5]);
  double a1 = atof(argv[6]), a2 = atof(argv[7]);
  unsigned r = atoi(argv[8]);
  ProdAdaKde2d kde(argv[3], h1, h2, a1, a2, r);
  kde.set_h1(h1); kde.set_h2(h2);
  double b1 = atof(argv[9]), e1 = atof(argv[10]), s1 = atof(argv[11]);
  double b2 = atof(argv[12]), e2 = atof(argv[13]), s2 = atof(argv[14]);

  // write 1d output
  ofstream f1out;
  open_for_writing(f1out, argv[2]);

  f1out << h1 << " " << h2 << " " << a1 << " " << a2 << endl;

  f1out << b1;
  for (double i = b1 + s1; i <= e1; i += s1) { f1out << " " << i; }
  f1out << endl;

  f1out << b2;
  for (double i = b2 + s2; i <= e2; i += s2) { f1out << " " << i; }
  f1out << endl;

  f1out << kde.f1(b1);
  for (double i = b1 + s1; i <= e1; i += s1) { f1out << " " << kde.f1(i); }
  f1out << endl;

  f1out << kde.f2(b2);
  for (double i = b2 + s2; i <= e2; i += s2) { f1out << " " << kde.f2(i); }
  f1out << endl;

  f1out.close();

  // write 2d output
  ofstream f2out;
  open_for_writing(f2out, argv[1]);

  f2out << b1;
  for (double i = b1 + s1; i <= e1; i += s1) { f2out << " " << i; }
  f2out << endl;

  f2out << b2;
  for (double i = b2 + s2; i <= e2; i += s2) { f2out << " " << i; }
  f2out << endl;

  for (double y = b2; y <= e2; y += s2) {
    f2out << kde(b1, y);
    for (double x = b1 + s1; x <= e1; x += s1) {
      f2out << " " << kde(x, y);
    }
    f2out << endl;
  }

  f2out.close();

  return 0;
}
