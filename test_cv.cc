#include <iostream>
#include <fstream>
#include <streambuf>

#include <file_io_utils.h>

#include "Kde2d.h"
#include "ProdKde2d.h"

using namespace std;

void cross_validate_scan(
    ProdKde2d &kde, 
    const vector<double> &candidates, 
    bool cv_x1,
    ostream &os) {

  double cv_min = 0.0, cv_argmin = 0.0;
  vector<double> results;
  for (auto &h : candidates) {
    kde.cv(results, h, cv_x1);
    if (results[2] < cv_min) { cv_min = results[2]; cv_argmin = h; }

    os << h << " ";
    os << results[1] << " ";
    os << results[2] << endl;
  }

  os << "best cv = " << cv_min << ", best h = " << cv_argmin << endl;
}

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


  string input_fname;
  vector<double> x1_candidates, x2_candidates;

  // evttype1
  input_fname = "data/evttype1.dat";
  x1_candidates = { 1e-1, 2e-2, 1e-2, 9e-3, 8e-3, 7.5e-3, 7e-3, 1e-3, 8e-4 };
  x2_candidates = { 1e-1, 9e-2, 8e-2, 7e-2, 6.5e-2, 6e-2, 5e-2, 3e-2 };
  ProdKde2d kde1(input_fname);

  os << "cross validating " << input_fname << ": "<< endl;
  os << "h1: " << endl;
  cross_validate_scan(kde1, x1_candidates, true, os);
  os << endl;
  os << "h2: " << endl;
  cross_validate_scan(kde1, x2_candidates, false, os);
  os << endl;

  // evttype2
  input_fname = "data/evttype2.dat";
  x1_candidates = { 1e-1, 2e-2, 1e-2, 9e-3, 8e-3, 7.5e-3, 7e-3, 1e-3, 8e-4 };
  x2_candidates = { 1e-1, 9e-2, 8e-2, 7e-2, 6.5e-2, 6e-2, 5e-2, 3e-2 };
  ProdKde2d kde2(input_fname);

  os << "cross validating " << input_fname << ": "<< endl;
  os << "h1: " << endl;
  cross_validate_scan(kde2, x1_candidates, true, os);
  os << endl;
  os << "h2: " << endl;
  cross_validate_scan(kde2, x2_candidates, false, os);
  os << endl;

  // evttype3
  input_fname = "data/evttype3_sub.dat";
  x1_candidates = { 7e-3, 6e-3, 5.5e-3, 5.2e-3, 5e-3 };
  x2_candidates = { 1e-1, 9e-2, 8e-2, 7e-2, 6e-2 };
  ProdKde2d kde3(input_fname);

  os << "cross validating " << input_fname << ": "<< endl;
  os << "h1: " << endl;
  cross_validate_scan(kde3, x1_candidates, true, os);
  os << endl;
  os << "h2: " << endl;
  cross_validate_scan(kde3, x2_candidates, false, os);
  os << endl;

  // evttype4
  input_fname = "data/evttype4_sub.dat";
  x1_candidates = { 6e-3, 5.7e-3, 5.5e-3, 5.2e-3 };
  x2_candidates = { 6e-2, 5.6e-2, 5.5e-2, 5e-2, 4e-2 };
  ProdKde2d kde4(input_fname);

  os << "cross validating " << input_fname << ": "<< endl;
  os << "h1: " << endl;
  cross_validate_scan(kde4, x1_candidates, true, os);
  os << endl;
  os << "h2: " << endl;
  cross_validate_scan(kde4, x2_candidates, false, os);
  os << endl;

  // evttype5
  input_fname = "data/evttype5_sub.dat";
  x1_candidates = { 4.35e-4, 4e-4, 3.5e-4, 2.5e-4, 1.5e-4 }; 
  x2_candidates = { 6e-2, 5e-2, 4e-2, 3.5e-2, 3e-2};
  ProdKde2d kde5(input_fname);

  os << "cross validating " << input_fname << ": "<< endl;
  os << "h1: " << endl;
  cross_validate_scan(kde5, x1_candidates, true, os);
  os << endl;
  os << "h2: " << endl;
  cross_validate_scan(kde5, x2_candidates, false, os);
  os << endl;

  if (fout.is_open()) {
    fout.close();
  }

  return 0;

}
