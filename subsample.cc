#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <random>

#include <file_io_utils.h>

using namespace std;

int main(int argc, char *argv[]) {

  if (argc < 3) {
    cerr << "usage: ./subsample input_fname " 
            "output_fname fraction" << endl;
    return 1;
  }

  double frac = atof(argv[3]);

  random_device rd;
  mt19937_64 e(rd());
  uniform_real_distribution<> d(0, 1);

  // read user input
  string line;
  ifstream fin; ofstream fout;
  open_for_reading(fin, argv[1]); open_for_writing(fout, argv[2]);
  while (getline(fin, line)) {
    if (d(e) < frac) {
      fout << line << "\n";
    }
  }

  return 0;
}
