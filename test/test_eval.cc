#include <iostream>
#include <fstream>
#include <streambuf>
#include <vector>
#include <utility>

#include <file_io_utils.h>

#include "Kde2d.h"
#include "ProdKde2d.h"

using namespace std;

int main(int argc, char *argv[]) {

  streambuf *buf = nullptr;
  ofstream fout;
  if (argc > 1) {
    open_for_writing(fout, argv[1]);
    buf = fout.rdbuf();
  } else {
    buf = cout.rdbuf();
  }
  ostream os(buf);

//opt_h1 = [ 1.97e-3, 1.75e-3, 3.92e-4, 5.75e-4, 1.09e-4 ]
//opt_h2 = [ 1.95e-2, 1.23e-2, 6.46e-3, 7.66e-3, 4.33e-3 ]

  // evttype1
  //string input_fname = "data/evttype1.dither.csv";
  //ProdKde2d kde(input_fname);
  //kde.set_h1(1.97e-3); kde.set_h2(1.95e-2);

  // evttype3
  string input_fname = "data/evttype3.dither.csv";
  ProdKde2d kde(input_fname);
  kde.set_h1(3.92e-4); kde.set_h2(6.46e-3);

  vector<pair<double, double>> results;
  kde.grid_evaluate_marginal(results, true, 20);
  for (auto e : results) {
    cout << e.first << " " << e.second << endl;
  }

  if (fout.is_open()) {
    fout.close();
  }

  return 0;

}
