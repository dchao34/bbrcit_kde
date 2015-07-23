#include <iostream>
#include <fstream>

#include <file_io_utils.h>

#include "ProdKde2d.h"

using namespace std;

int main() {

  ProdKde2d kde1("data/evttype1.dat");
  kde1.set_h1(0.008); kde1.set_h2(0.065);
  ProdKde2d kde2("data/evttype2.dat");
  kde2.set_h1(0.009); kde2.set_h2(0.065);
  ProdKde2d kde3("data/evttype3.dat");
  kde3.set_h1(0.0055); kde3.set_h2(0.07);
  ProdKde2d kde4("data/evttype4.dat");
  kde4.set_h1(0.0057); kde4.set_h2(0.055);
  ProdKde2d kde5("data/evttype5.dat");
  kde5.set_h1(0.000435); kde5.set_h2(0.04);

  ifstream fin; ofstream fout;
  open_for_reading(fin, "test_m.csv");
  open_for_writing(fout, "test_kdescore.csv");

  double x1, x2; int i = 0, c = 10000;
  while (fin >> x1 >> x2) {
    if (i % c == 0) {
      cout << "Processing lines " << i + 1 << " - " << i + c << endl;
    }

    fout << kde1(x1, x2) << " ";
    fout << kde2(x1, x2) << " ";
    fout << kde3(x1, x2) << " ";
    fout << kde4(x1, x2) << " ";
    fout << kde5(x1, x2);
    fout << endl;

    ++i;
  }

  fin.close(); fout.close();
  return 0;
}
