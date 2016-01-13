#include <iostream>
#include <vector>
#include <utility>
#include <random>

#include <Rectangle.h>
#include <Point.h>
#include <Kdtree.h>

using namespace std;
using bbrcit::Rectangle;
using bbrcit::Point;
using bbrcit::Kdtree;

int main() {

  random_device rd;
  default_random_engine g(rd());

  using Point2d = Point<2, double>;
  using Rectangle2d = Rectangle<2, double>;

  vector<Point2d> data;
  double start_x = -5, start_y = -5;
  double end_x = 5, end_y = 5;
  int N_x = 20, N_y = 20;
  double delta_x = (end_x - start_x) / N_x;
  double delta_y = (end_y - start_y) / N_y;
  for (int i = 0; i < N_x; ++i) {
    for (int j = 0; j < N_y; ++j) {
      data.push_back(Point2d({start_x + i*delta_x, start_y + j*delta_y}));
    }
  }

  for (const auto &p : data) { cout << p << " "; }
  cout << endl;
  cout << endl;
  //shuffle(data.begin(), data.end(), g);

  Kdtree<2,double> tree(data);

  Rectangle2d query({-1,-1}, {1,1});
  cout << query << endl;
  vector<Point2d> result;
  tree.range_search(result, query);
  for (const auto &p : result) { cout << p << " "; }
  cout << endl;


  /*using Point1d = Point<1, double>;
  using Rectangle1d = Rectangle<1, double>;

  vector<Point1d> data;
  double start_x = -5, end_x = 5;
  int N_x = 20;
  double delta_x = (end_x-start_x)/N_x;
  for (int i = 0; i < N_x; ++i) { data.push_back(Point1d({ start_x + i * delta_x })); }
  //shuffle(data.begin(), data.end(), g);

  Kdtree<1,double> tree(data);

  vector<Point1d> leaves; tree.report_leaves(leaves);
  for (auto l : leaves) { cout << l << " "; } cout << endl;

  Rectangle1d query({-2}, {2}); 
  cout << query << endl;
  vector<Point1d> result; tree.range_search(result, query);
  for (auto l : result) { cout << l << " "; } cout << endl;
  */

  return 0;
}
