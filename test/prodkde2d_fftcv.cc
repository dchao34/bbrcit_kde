#include <iostream>
#include <fstream>
#include <streambuf>
#include <cstdlib>

#include <file_io_utils.h>

#include <ProdKde2d.h>

using namespace std;

int main(int argc, char *argv[]) {

  if (argc < 3) {
    cerr << endl;
    cerr << "usage: ./prodkde2d_fftcv input_fname dim [output_fname]" << endl;
    cerr << "\tdim: 1 or 2. Cross validate first or second dimension. " << endl;
    cerr << endl;
    return EXIT_FAILURE;
  }

  streambuf *buf = nullptr;
  ofstream fout;
  if (argc > 3) {
    open_for_writing(fout, argv[3]);
    buf = fout.rdbuf();
  } else {
    buf = cout.rdbuf();
  }
  ostream os(buf);

  string input_fname(argv[1]);
  ProdKde2d kde(input_fname);

  bool cv_dim1 = (atoi(argv[2]) == 1) ?  true : false;

  vector<double> candidates;
  for (int i = 0; i < 200; i++) {
    candidates.push_back(0.5e-2 + i * 1e-3);
  }
  

  double cv_min = 0.0, cv_argmin = 0.0;
  vector<double> results;
  for (auto &h : candidates) {
    kde.fcv(results, h, 15, cv_dim1);
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
