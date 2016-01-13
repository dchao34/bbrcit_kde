#include <iostream>
#include <vector>
#include <utility>
#include <random>

#include <Point.h>
#include <Kdtree.h>

using namespace std;
using bbrcit::Point;
using bbrcit::Kdtree;

int main() {

  using Point2d = Point<2, double>;

  random_device rd;
  default_random_engine g(rd());

  vector<Point2d> data;
  double start_x = -5, start_y = -10;
  double end_x = 5, end_y = 10;
  int N_x = 100, N_y = 100;
  double delta_x = (end_x - start_x) / N_x;
  double delta_y = (end_y - start_y) / N_y;
  for (int i = 0; i < N_x; ++i) {
    for (int j = 0; j < N_y; ++j) {
      data.push_back(Point2d({start_x + i*delta_x, start_y + j*delta_y}));
    }
  }
  cout << data.size() << endl;
  shuffle(data.begin(), data.end(), g);
  Kdtree<2,double> tree(data);

  /*using Point1d = Point<1, double>;

  vector<Point1d> data;
  for (int i = 0; i < 100000; ++i) { data.push_back(Point1d({static_cast<double>(i)})); }
  shuffle(data.begin(), data.end(), g);
  for (int i = 0; i < 10; ++i) { cout << data[i] << " "; } cout << endl;
  cout << data.size() << endl;

  Kdtree<1,double> tree(data);
  vector<Point1d> leaves; tree.report_leaves(leaves);
  for (auto l : leaves) { cout << l << " "; } cout << endl;
  */

  return 0;
}
