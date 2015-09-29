#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <iostream>
#include <numeric>
#include <random>
#include <algorithm>
#include <cassert>

#include <file_io_utils.h>
#include <Kde2d.h>

using namespace std;

Kde2d::Kde2d(string data_fname, double bw1, double bw2) {

  h1 = bw1; h2 = bw2;

  ifstream fin;
  open_for_reading(fin, data_fname);

  // read the file line by line
  string line; double x1, x2;
  while (getline(fin, line)) {

    // read each line column by column
    istringstream sin(line);
    sin >> x1 >> x2;

    sample.push_back({x1, x2});
  }

  fin.close();
}

double Kde2d::operator()(double x1, double x2) {
  double result = 0.0;
  for (auto &p : sample) {
    result += gauss2d(x1, x2, p.first, p.second);
  }
  return result / sample.size();
}
