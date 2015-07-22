#include <iostream>
#include <fstream>

#include "Kde2d.h"
#include "ProdKde2d.h"

using namespace std;

int main() {

  ofstream fout("cv_results.txt");

  vector<double> h;

  // evttype1
  ProdKde2d kde1("data/evttype1.dat");
  fout << "evttype1: " << endl;
  fout << endl;
  h = { 1e-1, 2e-2, 1e-2, 9e-3, 8e-3, 7.5e-3, 7e-3, 1e-3, 8e-4 };
  kde1.cv(fout, h, true);
  fout << endl;
  h = { 1e-1, 9e-2, 8e-2, 7e-2, 6.5e-2, 6e-2, 5e-2, 3e-2 };
  kde1.cv(fout, h, false);
  fout << endl;

  // evttype2
  ProdKde2d kde2("data/evttype2.dat");
  fout << "evttype2: " << endl;
  fout << endl;
  h = { 1e-1, 2e-2, 1e-2, 9e-3, 8e-3, 7.5e-3, 7e-3, 1e-3, 8e-4 };
  kde2.cv(fout, h, true);
  fout << endl;
  h = { 1e-1, 9e-2, 8e-2, 7e-2, 6.5e-2, 6e-2, 5e-2, 3e-2 };
  kde2.cv(fout, h, false);
  fout << endl;

  // evttype3
  ProdKde2d kde3("data/evttype3_sub.dat");
  fout << "evttype3: " << endl;
  fout << endl;
  h = { 7e-3, 6e-3, 5.5e-3, 5.2e-3, 5e-3 };
  kde3.cv(fout, h, true);
  fout << endl;
  h = { 1e-1, 9e-2, 8e-2, 7e-2, 6e-2 };
  kde3.cv(fout, h, false);
  fout << endl;

  // evttype4
  ProdKde2d kde4("data/evttype4_sub.dat");
  fout << "evttype4: " << endl;
  fout << endl;
  h = { 6e-3, 5.7e-3, 5.5e-3, 5.2e-3 };
  kde4.cv(fout, h, true);
  fout << endl;
  h = { 6e-2, 5.6e-2, 5.5e-2, 5e-2, 4e-2 };
  kde4.cv(fout, h, false);
  fout << endl;

  // evttype5
  ProdKde2d kde5("data/evttype5_sub.dat");
  fout << "evttype5: " << endl;
  fout << endl;
  h = { 4.35e-4, 4e-4, 3.5e-4, 2.5e-4, 1.5e-4 }; 
  kde5.cv(fout, h, true);
  fout << endl;
  h = { 6e-2, 5e-2, 4e-2, 3.5e-2, 3e-2};
  kde5.cv(fout, h, false);
  fout << endl;

  return 0;

}
