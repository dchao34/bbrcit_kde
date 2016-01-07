#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <exception>
#include <utility>
#include <Point.h>
#include <Kdtree.h>

using namespace std;
using bbrcit::Point;
using bbrcit::Kdtree;

int main() {
  using Point1d = Point<1, double>;
  vector<Point1d> data;
  for (int i = 0; i < 10; ++i) { data.push_back(Point1d({static_cast<double>(i)})); }
  for (int i = 0; i < 10; ++i) { cout << data[i] << " "; } cout << endl;

  Kdtree<1,double> tree(data);
  tree.print();

  return 0;
}
