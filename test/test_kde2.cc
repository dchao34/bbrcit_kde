#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <utility>
#include <random>
#include <algorithm>

#include <KernelDensity.h>

#include "kde_test_utils.h"

using namespace std;
using bbrcit::KernelDensity;

using Kde1d = KernelDensity<1>;
using DataPoint1d = typename Kde1d::DataPointType;

int main() {

  int n_data = 10000;
  Kde1d kde;
  default_random_engine e;
  vector<DataPoint1d> data;

  int n_query = 100;
  vector<DataPoint1d> grid;
  generate_1dgrid(grid, -1, 1, n_query);

  double rel_err = 1e-6, abs_err = 1e-10;
  vector<DataPoint1d> direct_result, single_tree_result, dual_tree_result;

  ofstream fout;

  // test: unweighted
  data.clear();
  generate_bimodal_gaussian(e, data, n_data, 0.3, 0.1, -0.3, 0.1, 0.75);

  kde = Kde1d(data);
  kde.kernel().set_bandwidth(0.2);

  fout.open("test_kde2.csv");

  direct_result = grid;
  for (auto &q : direct_result) { kde.direct_eval(q); }
  std::sort(direct_result.begin(), direct_result.end(), 
      ReverseExactLexicoLess<DataPoint1d>);

  single_tree_result = grid;
  for (auto &q : single_tree_result) { kde.eval(q, rel_err, abs_err); }
  std::sort(single_tree_result.begin(), single_tree_result.end(), 
      ReverseExactLexicoLess<DataPoint1d>);

  dual_tree_result = grid;
  kde.eval(dual_tree_result, rel_err, abs_err);
  std::sort(dual_tree_result.begin(), dual_tree_result.end(), 
      ReverseExactLexicoLess<DataPoint1d>);

  for (int i = 0; i < n_query; ++i) {
    double direct = direct_result[i].attributes().value();
    double single = single_tree_result[i].attributes().value();
    double dual = dual_tree_result[i].attributes().value();

    double diff = std::abs(single - direct);
    if (diff > abs_err) {
      if (diff / std::max(single, direct) > rel_err) {
        cout << "relative precision lost for " << direct_result[i].point();
        cout << ": " << diff / std::max(single, direct); 
        cout << " expected " << rel_err << endl;
      } else {
        cout << "absolute precision lost for " << direct_result[i].point();
        cout << ": " << diff; 
        cout << " expected " << abs_err << endl;
      }
    }

    diff = std::abs(dual - direct);
    if (diff > abs_err) {
      if (diff / std::max(single, direct) > rel_err) {
        cout << "relative precision lost for " << direct_result[i].point();
        cout << ": " << diff / std::max(single, direct); 
        cout << " expected " << rel_err << endl;
      } else {
        cout << "absolute precision lost for " << direct_result[i].point();
        cout << ": " << diff; 
        cout << " expected " << abs_err << endl;
      }
    }

    fout << direct_result[i][0] << " ";
    fout << direct << " ";
    fout << single << " ";
    fout << dual << endl;
  }

  fout.close();

  // test: weighted. the proportions and weights dialed in for the 
  // bimodal gaussian below is such that the kde contains the "same"
  // probability content in either mode as `n_data` tends to infinity. 
  data.clear();
  generate_bimodal_gaussian(e, data, n_data, 0.3, 0.1, -0.3, 0.1, 0.25, 3);

  kde = Kde1d(data);
  kde.kernel().set_bandwidth(0.2);

  fout.open("test_kde2_weighted.csv");

  direct_result = grid;
  for (auto &q : direct_result) { kde.direct_eval(q); }
  std::sort(direct_result.begin(), direct_result.end(), 
      ReverseExactLexicoLess<DataPoint1d>);

  single_tree_result = grid;
  for (auto &q : single_tree_result) { kde.eval(q, rel_err, abs_err); }
  std::sort(single_tree_result.begin(), single_tree_result.end(), 
      ReverseExactLexicoLess<DataPoint1d>);

  dual_tree_result = grid;
  kde.eval(dual_tree_result, rel_err, abs_err);
  std::sort(dual_tree_result.begin(), dual_tree_result.end(), 
      ReverseExactLexicoLess<DataPoint1d>);

  for (int i = 0; i < n_query; ++i) {
    double direct = direct_result[i].attributes().value();
    double single = single_tree_result[i].attributes().value();
    double dual = dual_tree_result[i].attributes().value();

    double diff = std::abs(single - direct);
    if (diff > abs_err) {
      if (diff / std::max(single, direct) > rel_err) {
        cout << "relative precision lost for " << direct_result[i].point();
        cout << ": " << diff / std::max(single, direct); 
        cout << " expected " << rel_err << endl;
      } else {
        cout << "absolute precision lost for " << direct_result[i].point();
        cout << ": " << diff; 
        cout << " expected " << abs_err << endl;
      }
    }

    diff = std::abs(dual - direct);
    if (diff > abs_err) {
      if (diff / std::max(single, direct) > rel_err) {
        cout << "relative precision lost for " << direct_result[i].point();
        cout << ": " << diff / std::max(single, direct); 
        cout << " expected " << rel_err << endl;
      } else {
        cout << "absolute precision lost for " << direct_result[i].point();
        cout << ": " << diff; 
        cout << " expected " << abs_err << endl;
      }
    }

    fout << direct_result[i][0] << " ";
    fout << direct << " ";
    fout << single << " ";
    fout << dual << endl;
  }

  fout.close();

  return 0;
}
