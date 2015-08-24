#include <iostream>
#include <fstream>
#include <streambuf>

#include <file_io_utils.h>

#include "ProdAdaKde2d.h"

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

  string input_fname;
  vector<double> h1_candidates;
  vector<double> h2_candidates;
  vector<double> results;
  double cv_min, cv_argmin;

  // evttype1
  input_fname = "data/evttype1.cv.csv";

  ProdAdaKde2d kde(input_fname, 1.97e-3, 1.95e-2, 0.5, 0.5, 20);

  h1_candidates = { 3e-3, 2e-3, 1.5e-3, 1e-3, 9e-4 };
  h2_candidates = { 3e-2, 2e-2, 1.5e-2, 1e-2, 9e-3 };

  cout << "h1: " << endl;
  cv_min = 0.0, cv_argmin = 0.0;
  for (auto &h : h1_candidates) {
    kde.cv(results, h, true);
    if (results[3] < cv_min) { cv_min = results[3]; cv_argmin = h; }

    os << h << " ";
    os << results[1] << " ";
    os << results[2] << " ";
    os << results[3] << endl;
  }

  os << "best cv = " << cv_min << ", best h = " << cv_argmin << endl;
  cout << endl;

  cout << "h2: " << endl;
  cv_min = 0.0, cv_argmin = 0.0;
  for (auto &h : h2_candidates) {
    kde.cv(results, h, false);
    if (results[3] < cv_min) { cv_min = results[3]; cv_argmin = h; }

    os << h << " ";
    os << results[1] << " ";
    os << results[2] << " ";
    os << results[3] << endl;
  }

  os << "best cv = " << cv_min << ", best h = " << cv_argmin << endl;
  cout << endl;


  // evttype2
  /*input_fname = "data/evttype2.cv.csv";

  ProdAdaKde2d kde(input_fname, 1.75e-3, 1.23e-2, 0.5, 0.5, 20);

  h1_candidates = { 1e-1, 1e-2, 1e-3, 1e-4 };
  h2_candidates = { 1e-1, 1e-2, 1e-3, 1e-4 };

  cout << "h1: " << endl;
  cv_min = 0.0, cv_argmin = 0.0;
  for (auto &h : h1_candidates) {
    kde.cv(results, h, true);
    if (results[3] < cv_min) { cv_min = results[3]; cv_argmin = h; }

    os << h << " ";
    os << results[1] << " ";
    os << results[2] << " ";
    os << results[3] << endl;
  }

  os << "best cv = " << cv_min << ", best h = " << cv_argmin << endl;
  cout << endl;

  cout << "h2: " << endl;
  cv_min = 0.0, cv_argmin = 0.0;
  for (auto &h : h2_candidates) {
    kde.cv(results, h, false);
    if (results[3] < cv_min) { cv_min = results[3]; cv_argmin = h; }

    os << h << " ";
    os << results[1] << " ";
    os << results[2] << " ";
    os << results[3] << endl;
  }

  os << "best cv = " << cv_min << ", best h = " << cv_argmin << endl;
  cout << endl;
  */
  

  if (fout.is_open()) {
    fout.close();
  }

  return 0;

}
