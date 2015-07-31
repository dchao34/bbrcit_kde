#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <random>

#include "file_io_utils.h"

using namespace std;

int main(int argc, char *argv[]) {

  if (argc < 3) {
    cerr << "usage: ./dither input_fname output_fname" << endl;
    return 1;
  }

  random_device rd;
  mt19937_64 e(rd());
  uniform_real_distribution<> d(-0.5, 0.5);

  string line;
  ifstream fin; ofstream fout;
  open_for_reading(fin, argv[1]); open_for_writing(fout, argv[2]);
  while (getline(fin, line)) {
    stringstream ss(line);
    double c0, c1; ss >> c0; ss >> c1;
    c0 += d(e) * 1e-3; c1 += d(e) * 1e-3;
    fout << c0 << " " << c1 << endl;
  }

  return 0;
}
