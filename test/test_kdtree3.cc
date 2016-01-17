#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <random>

#include <Kdtree.h>

using namespace std;
using bbrcit::Kdtree;

using Kdtree1d = Kdtree<1>;
using Point1d = typename Kdtree1d::PointType;
using Rectangle1d = typename Kdtree1d::RectangleType;

int main() {

  vector<Point1d> points1;
  for (int i = 0; i < 10; ++i) { points1.push_back({{static_cast<double>(i)}}); }
  for (int i = 0; i < 10; ++i) { points1.push_back({{static_cast<double>(i)}}); }
  for (const auto &p : points1) { cout << p << endl; }
  cout << endl;

  Kdtree1d tr1(points1);
  tr1.report_leaves(cout);

  cout << endl;

  cout << points1.size() << " -> " << tr1.size() << endl;

  cout << endl;


  return 0;
}
