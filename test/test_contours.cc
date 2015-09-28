#include <iostream>
#include <fstream>
#include <sstream>

#include <Kde2d.h>
#include <ProdKde2d.h>

using namespace std;

int main() {

  double h1 = 0, h2 = 0;
  string fname; ofstream fout;

  // evttype1
  cout << "evttype1" << endl;
  ProdKde2d kde1("evttype1.dat");
  h1 = 0.008, h2 = 0.065;
  fname = "contour/evttype1_cv.dat";

  kde1.set_h1(h1); kde1.set_h2(h2);
  fout.open(fname);
  for (double x1 = 0; x1 < 0.1; x1 += 0.001) {
    for (double x2 = 0; x2 < 1.0; x2 += 0.01) {
      fout << x1 << " " << x2 << " " << kde1(x1, x2) << "\n";
    }
  }
  fout.close();
  cout << "h1 = " << h1;
  cout << ", h2 = " << h2;
  cout << ", saved to " << fname << endl;


  // evttype2
  cout << "evttype2" << endl;
  ProdKde2d kde2("evttype2.dat");
  h1 = 0.009, h2 = 0.065;
  fname = "contour/evttype2_cv.dat";

  kde2.set_h1(h1); kde2.set_h2(h2);
  fout.open(fname);
  for (double x1 = 0; x1 < 0.1; x1 += 0.001) {
    for (double x2 = 0; x2 < 1.0; x2 += 0.01) {
      fout << x1 << " " << x2 << " " << kde2(x1, x2) << "\n";
    }
  }
  fout.close();
  cout << "h1 = " << h1;
  cout << ", h2 = " << h2;
  cout << ", saved to " << fname << endl;


  // evttype3
  cout << "evttype3" << endl;
  ProdKde2d kde3("evttype3.dat");
  h1 = 0.0055, h2 = 0.07;
  fname = "contour/evttype3_cv.dat";

  kde3.set_h1(h1); kde3.set_h2(h2);
  fout.open(fname);
  for (double x1 = 0; x1 < 0.1; x1 += 0.001) {
    for (double x2 = 0; x2 < 1.0; x2 += 0.01) {
      fout << x1 << " " << x2 << " " << kde3(x1, x2) << "\n";
    }
  }
  fout.close();
  cout << "h1 = " << h1;
  cout << ", h2 = " << h2;
  cout << ", saved to " << fname << endl;


  // evttype4
  cout << "evttype4" << endl;
  ProdKde2d kde4("evttype4.dat");
  h1 = 0.0057, h2 = 0.055;
  fname = "contour/evttype4_cv.dat";

  kde4.set_h1(h1); kde4.set_h2(h2);
  fout.open(fname);
  for (double x1 = 0; x1 < 0.1; x1 += 0.001) {
    for (double x2 = 0; x2 < 1.0; x2 += 0.01) {
      fout << x1 << " " << x2 << " " << kde4(x1, x2) << "\n";
    }
  }
  fout.close();
  cout << "h1 = " << h1;
  cout << ", h2 = " << h2;
  cout << ", saved to " << fname << endl;

  // evttype5
  cout << "evttype5" << endl;
  ProdKde2d kde5("evttype5.dat");
  h1 = 0.000435, h2 = 0.04;
  fname = "contour/evttype5_cv.dat";

  kde5.set_h1(h1); kde5.set_h2(h2);
  fout.open(fname);
  for (double x1 = 0; x1 < 0.1; x1 += 0.001) {
    for (double x2 = 0; x2 < 1.0; x2 += 0.01) {
      fout << x1 << " " << x2 << " " << kde5(x1, x2) << "\n";
    }
  }
  fout.close();
  cout << "h1 = " << h1;
  cout << ", h2 = " << h2;
  cout << ", saved to " << fname << endl;

  return 0;
}
