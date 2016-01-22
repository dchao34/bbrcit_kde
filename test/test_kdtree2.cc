#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <random>

#include <Kdtree.h>

using namespace std;
using bbrcit::Kdtree;

using Kdtree1d = Kdtree<1>;
using DataPoint1d = typename Kdtree1d::DataPointType;
using Rectangle1d = typename Kdtree1d::RectangleType;

int main() {
  cout << endl;

  // test: kdtree construction using different leaf_nmax
  cout << "+ kdtree construction with different leaf_max. " << endl;
  cout << endl;
  vector<DataPoint1d> points1d;
  for (int i = 0; i < 8; ++i) { points1d.push_back({{static_cast<double>(i)}}); }
  Kdtree1d tr1(points1d, 2);
  vector<pair<size_t, size_t>> leaves;
  tr1.report_leaves(leaves);
  cout << "  output:  ";
  for (const auto &l : leaves) { cout << "[" << l.first << ", " << l.second << "] "; }
  cout << endl;
  cout << "  compare: [0, 1] [2, 3] [4, 5] [6, 7]" << endl;
  cout << endl;

  Kdtree1d tr2(points1d, 1);
  leaves.clear(); tr2.report_leaves(leaves);
  cout << "  output:  ";
  for (const auto &l : leaves) { cout << "[" << l.first << ", " << l.second << "] "; }
  cout << endl;
  cout << "  compare: [0, 0] [1, 1] [2, 2] [3, 3] [4, 4] [5, 5] [6, 6] [7, 7]" << endl;
  cout << endl;

  Kdtree1d tr3(points1d, 5);
  leaves.clear(); tr3.report_leaves(leaves);
  cout << "  output:  ";
  for (const auto &l : leaves) { cout << "[" << l.first << ", " << l.second << "] "; }
  cout << endl;
  cout << "  compare: [0, 3] [4, 7]" << endl;
  
  cout << endl;

  // test: duplicate keys
  cout << "+ kdtree with duplicate keys. " << endl;
  cout << endl;

  cout << "  output:";
  points1d.clear();
  for (int i = 0; i < 4; ++i) { points1d.push_back({{static_cast<double>(i)}}); }
  for (int i = 0; i < 4; ++i) { points1d.push_back({{static_cast<double>(i)}}); }
  Kdtree1d tr4(points1d, 2);
  for (const auto &p : tr4.points()) { 
    cout << "  " << p;
  }
  cout << endl;

  cout << "  compare: { (0), (2) }  { (1), (2) }  { (2), (2) }  { (3), (2) }";
  cout << endl;

  cout << endl;


  return 0;
}
