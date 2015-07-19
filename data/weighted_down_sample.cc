#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>

#include "file_io_utils.h"

using namespace std;

// --------------------------------------------------------------------
// Input: 
//     File with one record per line, delimited with ",". One of the
//     column entry specifies a record weight.
// Output: 
//     File with one record per line, delimited with ",". It is a 
//     subset of the input file, down sampled according to record 
//     weight.
// --------------------------------------------------------------------

int main(int argc, char *argv[]) {

  if (argc != 3) {
    cerr << "usage: ./weighted_down_sample original_fname down_sampled_fname" << endl;
    return 1;
  }

  // open files for reading/writing
  ifstream fin; ofstream fout; 
  open_for_reading(fin, argv[1]);
  open_for_writing(fout, argv[2]);

  // write title line
  string line; getline(fin, line); 
  fout << line << "\n";

  // down sample according to weight
  random_device rd;
  mt19937_64 e(rd());
  uniform_real_distribution<> u(0, 1);
  size_t line_cnt = 0;
  vector<string> tokens;
  while (getline(fin, line)) {
    tokenize(line, tokens, ",");
    if (u(e) <= stod(tokens[3])) {
      fout << line << "\n";
    }
    ++line_cnt;
  }
  
  // clean up
  fin.close(); fout.close();
  cout << "Read " << line_cnt << " lines." << endl;

  return 0;
}
