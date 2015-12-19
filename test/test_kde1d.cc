#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <Point.h>
#include <Kde1d.h>

using namespace std;

int main() {

  // read generated data into a vector
  string input_fname("bimodal_gauss1d.csv"); ifstream fin(input_fname);
  cout << "reading generated data from " << input_fname << endl;
  vector<Point1d> sample;
  double value; while (fin >> value) { sample.push_back(value); }

  // evaluate kde at various bandwidths 
  Kde1d kde(sample); int n = 2000; double bound = 2.0;

  cout << "evaluating kde with bandwidth = " << kde.get_bandwidth() << endl;
  vector<double> result_bw1;
  for (int i = 0; i < n; ++i) {
    int x = -bound + 2*bound*i/n;
    result_bw1.push_back(kde.eval(x));
  }

  vector<double> result_bw2;
  kde.set_bandwidth(0.01);
  cout << "evaluating kde with bandwidth = " << kde.get_bandwidth() << endl;
  for (int i = 0; i < n; ++i) {
    int x = -bound + 2*bound*i/n;
    result_bw2.push_back(kde.eval(x));
  }

  vector<double> result_bw3;
  kde.set_bandwidth(1e-7);
  cout << "evaluating kde with bandwidth = " << kde.get_bandwidth() << endl;
  for (int i = 0; i < n; ++i) {
    int x = -bound + 2*bound*i/n;
    result_bw3.push_back(kde.eval(x));
  }

  
  // evaluate kde at various bandwidths 
  string output_fname("kde1d_gauss1d.csv"); ofstream fout(output_fname);
  cout << "saving results to " << output_fname << endl;
  for (int i = 0; i < n; ++i) {
    double x = -bound + 2*bound*i/n;
    fout << x << " ";
    fout << result_bw1[i] << " ";
    fout << result_bw2[i] << " ";
    fout << result_bw3[i] << endl;
  }

  cout << "done. "<< endl;

  return 0;
}
