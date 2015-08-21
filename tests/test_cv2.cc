#include <iostream>
#include <fstream>
#include <streambuf>

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

  // evttype1
  string input_fname = "data/evttype1.cv.csv";
  ProdKde2d kde(input_fname);

  vector<double> candidates;
  for (int i = 0; i < 10; i++) {
    candidates.push_back(1e-3 + i * 1e-3);
  }
  

  double cv_min = 0.0, cv_argmin = 0.0;
  vector<double> results;
  for (auto &h : candidates) {
    kde.cv(results, h, true);
    if (results[2] < cv_min) { cv_min = results[2]; cv_argmin = h; }

    os << h << " ";
    os << results[1] << " ";
    os << results[2] << endl;
  }

  os << "best cv = " << cv_min << ", best h = " << cv_argmin << endl;


  if (fout.is_open()) {
    fout.close();
  }

  return 0;

}
