#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <string>

#include <file_io_utils.h>

#include "ProdAdaKde2d.h"

using namespace std;

using vp = vector<pair<double,double>>;
using sample_no = vp::size_type;

int main(int argc, char *argv[]) {

  if (argc < 5) {
    cerr << "usage: ./cache_akde_score datadir output_file_name first_idx last_idx" << endl;
    cerr << endl;
    cerr << "first_idx and last_idx are 0 indexed." << endl;
    cerr << "processes the inclusive range [first_idx, last_idx]." << endl;
    return 1;
  }

  ProdAdaKde2d kde1(string(argv[1]) + "/evttype1.dither.csv", 1.97e-3, 1.95e-2, 0.28125, 0.28125, 20);
  ProdAdaKde2d kde2(string(argv[1]) + "/evttype2.dither.csv", 1.75e-3, 1.23e-2, 0.28125, 0.28125, 20);
  ProdAdaKde2d kde3(string(argv[1]) + "/evttype3.dither.csv", 3.92e-4, 6.46e-3, 0.28125, 0.28125, 20);
  ProdAdaKde2d kde4(string(argv[1]) + "/evttype4.dither.csv", 5.75e-4, 7.66e-3, 0.28125, 0.28125, 20);
  ProdAdaKde2d kde5(string(argv[1]) + "/evttype5.dither.csv", 1.09e-4, 4.33e-3, 0.28125, 0.28125, 20);

  ifstream fin; open_for_reading(fin, string(argv[1]) + "/test.dither.csv");
  vp points;
  double x1, x2;
  while (fin >> x1 >> x2) {
    points.push_back(make_pair(x1, x2));
  }
  fin.close();

  ofstream fout;
  open_for_writing(fout, argv[2]);
  sample_no beg_idx = atoi(argv[3]), end_idx = atoi(argv[4]) + 1;
  sample_no print_frequency = 5;
  for (auto i = beg_idx; i < end_idx; ++i) {
    if (i % print_frequency == 0) {
      cout << "Processing lines " << i + 1 << " - " << i + print_frequency << endl;
    }
    fout << kde1(points[i].first, points[i].second) << " ";
    fout << kde2(points[i].first, points[i].second) << " ";
    fout << kde3(points[i].first, points[i].second) << " ";
    fout << kde4(points[i].first, points[i].second) << " ";
    fout << kde5(points[i].first, points[i].second);
    fout << endl;
  }
  fout.close();

  return 0;
}
