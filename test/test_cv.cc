#include <iostream>
#include <fstream>
#include <streambuf>

#include <file_io_utils.h>

#include <Kde2d.h>
#include <ProdKde2d.h>

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
  input_fname = "data/evttype1.cv.csv";
  x1_candidates = { 2.35e-3, 2.31e-3, 2.29e-3, 2.27e-3, 2.22e-3};
  x2_candidates = { 2.19e-3, 2.18e-2, 2.16e-2, 2.15e-2, 2.1e-3 };
  ProdKde2d kde1(input_fname);

  os << "cross validating " << input_fname << ": "<< endl;
  os << "h1: " << endl;
  cross_validate_scan(kde1, x1_candidates, true, os);
  os << endl;
  os << "h2: " << endl;
  cross_validate_scan(kde1, x2_candidates, false, os);
  os << endl;

  // evttype2
  input_fname = "data/evttype2.cv.csv";
  x1_candidates = { 2.14e-3, 2.12e-3, 2.1e-3, 2.05e-3, 2.0e-3 };
  x2_candidates = { 1.7e-2, 1.66e-2, 1.64e-2, 1.63e-2, 1.6e-2 };
  ProdKde2d kde2(input_fname);

  os << "cross validating " << input_fname << ": "<< endl;
  os << "h1: " << endl;
  cross_validate_scan(kde2, x1_candidates, true, os);
  os << endl;
  os << "h2: " << endl;
  cross_validate_scan(kde2, x2_candidates, false, os);
  os << endl;


  // evttype3
  input_fname = "data/evttype3.cv.csv";
  x1_candidates = { 1.3e-3, 1.26e-3, 1.24e-3, 1.22e-3, 1.20e-3 };
  x2_candidates = { 2.2e-2, 2.0e-2, 1.91e-2, 1.8e-2, 1.6e-2 };
  ProdKde2d kde3(input_fname);

  os << "cross validating " << input_fname << ": "<< endl;
  os << "h1: " << endl;
  cross_validate_scan(kde3, x1_candidates, true, os);
  os << endl;
  os << "h2: " << endl;
  cross_validate_scan(kde3, x2_candidates, false, os);
  os << endl;

  // evttype4
  input_fname = "data/evttype4.cv.csv";
  x1_candidates = { 1.47e-3, 1.44e-3, 1.42e-3, 1.41e-3, 1.39e-3};
  x2_candidates = { 1.7e-2, 1.64e-2, 1.63e-2, 1.62e-2, 1.58e-2 };
  ProdKde2d kde4(input_fname);

  os << "cross validating " << input_fname << ": "<< endl;
  os << "h1: " << endl;
  cross_validate_scan(kde4, x1_candidates, true, os);
  os << endl;
  os << "h2: " << endl;
  cross_validate_scan(kde4, x2_candidates, false, os);
  os << endl;

  // evttype5
  input_fname = "data/evttype5.cv.csv";
  x1_candidates = { 1.69e-4, 1.66e-4, 1.65e-4, 1.64e-4, 1.62e-4}; 
  x2_candidates = { 1.5e-2, 1.46e-2, 1.44e-2, 1.42e-2, 1.4e-2};
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
