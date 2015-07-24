#include <iostream>
#include <iomanip>

#include <file_io_utils.h>

#include "ProdKde2d.h"

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
//     *f1* is the marginal kde for the 1st dimension.
//     *f2* is the marginal kde for the 2nd dimension.
//     *m* is the total number of grid points over the 1st dimension.
//     *n* is the total number of grid points over the 2nd dimension.
//     line 1: b1 b1+s1 b1+2*s1 ... b1+(m-1)*s1
//     line 2: b2 b2+s2 b1+2*s2 ... b1+(n-1)*s2
//     line 3: f1(b1) f1(b1+s1) ... f1(b1+(n-1)*s1)
//     line 4: f2(b2) f2(b2+s2) ... f2(b2+(n-2)*s2)
int main(int argc, char *argv[]) {

  if (argc < 7) {
    cerr << endl;
    cerr << "usage: ./kde_scan 2d_out 1d_out data_in bw1 bw2 b1 e1 s1 b2 e2 s2" << endl;
    cerr << endl;
    cerr << left << setfill('.');
    cerr << setw(30) << "1. 2d_out" <<  "file to save 2d scan. " << endl;
    cerr << setw(30) << "2. 1d_out" << "file to save 1d scans. " << endl;
    cerr << setw(30) << "3. data_in" << "file containing the sample points. " << endl;
    cerr << setw(30) << "4-5. bw1 bw2" << "bandwidths in 1st and 2nd dimensions. " << endl;
    cerr << setw(30) << "6-8. b1, e1, s1" << "grid points in 1st dimension. " << endl;
    cerr << setw(30) << setfill(' ') << "" << "covers the points [<b1>, <e1>) "
                                              "in *increments* of <s1>." << endl; 
    cerr << setfill('.');
    cerr << setw(30) << "9-11. b2, e2, s2" << "grid points in 2nd dimension. " << endl;
    cerr << endl;
    return 1;
  }

  // load kde
  ProdKde2d kde(argv[3]);
  kde.set_h1(atof(argv[4])); kde.set_h2(atof(argv[5]));
  double b1 = atof(argv[6]), e1 = atof(argv[7]), s1 = atof(argv[8]);
  double b2 = atof(argv[9]), e2 = atof(argv[10]), s2 = atof(argv[11]);

  // write 1d output
  ofstream f1out;
  open_for_writing(f1out, argv[2]);

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
