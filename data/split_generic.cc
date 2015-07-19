#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <random>

#include "file_io_utils.h"

using namespace std;

// --------------------------------------------------------------------
// Splits input samples randomly between training and testing files. 
// --------------------------------------------------------------------

int main() {

  // open files for reading/writing
  ifstream fin; ofstream tr_fout, te_fout;
  open_for_reading(fin, "fitdata_generic.csv");
  open_for_writing(tr_fout, "train_generic.csv");
  open_for_writing(te_fout, "test_generic.csv");

  // write title line
  string line; getline(fin, line); 
  tr_fout << line << "\n"; te_fout << line << "\n";


  // split lines between the two files
  random_device rd;
  mt19937_64 e(rd());
  uniform_real_distribution<> u(0, 1);
  size_t line_cnt = 0;
  while (getline(fin, line)) {
    ++line_cnt;
    if (u(e) > 0.714) {
      tr_fout << line << "\n";
    } else {
      te_fout << line << "\n";
    }
  }
  
  fin.close(); tr_fout.close(); te_fout.close();
  cout << "Read/divided " << line_cnt << " lines." << endl;

  return 0;
}
