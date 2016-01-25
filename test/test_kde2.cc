#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <random>

#include <KernelDensity.h>

using namespace std;
using bbrcit::KernelDensity;

using Kde1d = KernelDensity<1>;
using DataPoint1d = typename Kde1d::DataPointType;
using Rectangle1d = typename Kde1d::GeomRectangleType;

int main() {
  cout << endl;

  // test: kde construction using different leaf_nmax
  cout << "+ kde construction with different leaf_max. " << endl;
  cout << endl;
  vector<DataPoint1d> points1d;
  for (int i = 0; i < 8; ++i) { points1d.push_back({{static_cast<double>(i)}}); }
  Kde1d tr1(points1d, 1, 2);
  vector<pair<size_t, size_t>> leaves;
  tr1.report_leaves(leaves);
  cout << "  output:  ";
  for (const auto &l : leaves) { cout << "[" << l.first << ", " << l.second << "] "; }
  cout << endl;
  cout << "  compare: [0, 1] [2, 3] [4, 5] [6, 7]" << endl;
  cout << endl;

  Kde1d tr2(points1d, 1, 1);
  leaves.clear(); tr2.report_leaves(leaves);
  cout << "  output:  ";
  for (const auto &l : leaves) { cout << "[" << l.first << ", " << l.second << "] "; }
  cout << endl;
  cout << "  compare: [0, 0] [1, 1] [2, 2] [3, 3] [4, 4] [5, 5] [6, 6] [7, 7]" << endl;
  cout << endl;

  Kde1d tr3(points1d, 1, 5);
  leaves.clear(); tr3.report_leaves(leaves);
  cout << "  output:  ";
  for (const auto &l : leaves) { cout << "[" << l.first << ", " << l.second << "] "; }
  cout << endl;
  cout << "  compare: [0, 3] [4, 7]" << endl;
  
  cout << endl;

  // test: duplicate keys
  cout << "+ kde with duplicate keys. " << endl;
  cout << endl;

  cout << "  output:";
  points1d.clear();
  for (int i = 0; i < 4; ++i) { points1d.push_back({{static_cast<double>(i)}}); }
  for (int i = 0; i < 4; ++i) { points1d.push_back({{static_cast<double>(i)}}); }
  Kde1d tr4(points1d, 1);
  for (const auto &p : tr4.data_points()) { 
    cout << "  " << p;
  }
  cout << endl;

  cout << "  compare: { (0), (2) }  { (1), (2) }  { (2), (2) }  { (3), (2) }";
  cout << endl;

  cout << endl;

  // test: root attributes
  points1d.clear();
  for (int i = 0; i < 4; ++i) { points1d.push_back({{1.0*i}, {0.1*i}}); }
  Kde1d tr5(points1d, 1);
  cout << "+ kde root attributes. " << endl;
  cout << endl;
  cout << "  output:  " << tr5.root_attributes() << " (c.f. (0.6))" << endl;
  cout << endl;

  // test: range search
  points1d.clear();
  for (int i = 0; i < 1000; ++i) { points1d.push_back({{0.001*i}}); }
  Kde1d tr6(points1d, 1);
  vector<DataPoint1d> result;
  Rectangle1d query({0.990}, {0.995});
  tr6.range_search(query, result);
  cout << "+ kde range search. " << endl;
  cout << endl;
  cout << "  output:  "; for (const auto &p : result) { cout << " " << p[0]; } cout << ", " << tr6.root_attributes() << endl; 
  cout << "  compare:  0.99 0.991 0.992 0.993 0.994 0.995, (1000)" << endl;
  cout << endl;


  return 0;
}
