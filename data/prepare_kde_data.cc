#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>

#include "file_io_utils.h"

using namespace std;

// -------------------------------------------------------------------------
// Input:
//   1. Input csv file name. 
//
//      - First row is the title. It must be in the following format:
//
//        rf_useopt_score,rf_dvsdstar_sigmc_score,tag_lp3,event_weight,mc_evttype
//
//      - Subsequent rows are individual data records, whose columns are 
//        ordered to correspond to the name specified in the title. 
//
//   2. Output data file name. 
//
//   3. 2 positive integers. Column indices for the two features. 
//
//      It is your responsibility to give valid indices. 
//
//   4. At least 1 positive integer. These mc_evttype's will be 
//      extracted to the output file.
//
//      It is your responsibility to give valid component indices. 
// -------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  if (argc < 6) {
    cerr << "usage: ./prepare_kde_data input_fname "
            "output_fname col1 col2 mc_evttype ... " << endl;
    return 1;
  }

  // read user input
  ifstream fin; ofstream fout;
  open_for_reading(fin, argv[1]); open_for_writing(fout, argv[2]);
  int col1 = atoi(argv[3]), col2 = atoi(argv[4]);
  vector<int> evttypes;
  for (int i = 5; i < argc; i++) {
    evttypes.push_back(atoi(argv[i]));
  }

  // stream records in a line at a time.
  string line; getline(fin, line);
  vector<string> columns;
  while (getline(fin, line)) {
    tokenize(line, columns, ",");
    
    // check whether the line has the correct type and save the specified columns.
    if (find(evttypes.begin(), evttypes.end(), 
          stoi(columns.back())) != evttypes.end()) {
      fout << columns[col1] << " " << columns[col2] << "\n";
    }
  }

  fin.close(); fout.close();
  return 0;
}
