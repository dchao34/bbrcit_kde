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

using Point2d = Point<2, double>;
using Rectangle2d = Rectangle<2, double>;

vector<Point2d> generate_grid(
    double start_x, double end_x, 
    double start_y, double end_y, 
    int N_x, int N_y) {

  vector<Point2d> result;

  double delta_x = (end_x - start_x) / N_x, delta_y = (end_y - start_y) / N_y;
  for (int i = 0; i < N_x; ++i) { 
    for (int j = 0; j < N_y; ++j) {
      result.push_back(Point2d({start_x + i*delta_x, start_y + j*delta_y}));
    }
  }

  return result;
}

int main() {

  cout << endl;

  // test: default constructor
  Kdtree<2,double> tr0;
  cout << "+ Default constructor: " << tr0.empty() << " (cf. 1)" << endl;
  cout << endl;

  // test: vector constructor
  vector<Point2d> points1 = generate_grid(0, 1, 0, 1, 2, 2);
  Kdtree<2, double> tr1(points1);
  vector<Point2d> leaves; tr1.report_leaves(leaves);
  cout << "+ vector constructor (1): { ";
  for (const auto &p : leaves) { cout << p << " "; }
  cout << "} (c.f. { (0, 0) (0, 0.5) (0.5, 0) (0.5, 0.5) })" << endl;

  vector<Point2d> points2 = generate_grid(0, -1, 0, -1, 2, 2);
  Kdtree<2, double> tr2(std::move(points2));
  leaves.clear(); tr2.report_leaves(leaves);
  cout << "+ vector constructor (2): { ";
  for (const auto &p : leaves) { cout << p << " "; }
  cout << "} (c.f. { (-0.5, -0.5) (-0.5, 0) (0, -0.5), (0, 0) })" << endl;

  cout << endl;

  // test: copy constructor
  Kdtree<2,double> *p_tr = new Kdtree<2,double>(generate_grid(0, 1, 0, -1, 2, 2));
  Kdtree<2,double> tr3(*p_tr);
  delete p_tr;
  leaves.clear(); tr3.report_leaves(leaves);
  cout << "+ copy constructor (1): { ";
  for (const auto &p : leaves) { cout << p << " "; }
  cout << "} (c.f. { (0, -0.5) (0, 0) (0.5, -0.5) (0.5, 0) })" << endl;

  Kdtree<2,double> temp_tr(generate_grid(0, -1, 0, 1, 2, 2));
  Kdtree<2,double> tr4(std::move(temp_tr));
  leaves.clear(); tr3.report_leaves(leaves);
  cout << "+ move constructor (1): { ";
  for (const auto &p : leaves) { cout << p << " "; } cout << "}, ";
  cout << temp_tr.empty(); 
  cout << " (c.f. { (-0.5, 0) (-0.5, -0.5) (0, 0) (0, 0.5) }, 1)" << endl;
  cout << endl;

  // test: swap
  cout << "+ swap (1): " << endl;;
  swap(tr3, tr4); 
  leaves.clear(); tr3.report_leaves(leaves);
  cout << "\t { ";
  for (const auto &p : leaves) { cout << p << " "; } cout << "}, ";
  cout << "(c.f. { (-0.5, 0) (-0.5, -0.5) (0, 0) (0, 0.5) })" << endl;
  leaves.clear(); tr4.report_leaves(leaves);
  cout << "\t { ";
  for (const auto &p : leaves) { cout << p << " "; } cout << "}, ";
  cout << "(c.f. { (0, -0.5) (0, 0) (0.5, -0.5) (0.5, 0) })" << endl;
  swap(tr3, tr4); 
  leaves.clear(); tr3.report_leaves(leaves);
  cout << "\t { ";
  for (const auto &p : leaves) { cout << p << " "; } cout << "}, ";
  cout << "(c.f. { (0, -0.5) (0, 0) (0.5, -0.5) (0.5, 0) })" << endl;
  leaves.clear(); tr4.report_leaves(leaves);
  cout << "\t { ";
  for (const auto &p : leaves) { cout << p << " "; } cout << "}, ";
  cout << "(c.f. { (-0.5, 0) (-0.5, -0.5) (0, 0) (0, 0.5) })" << endl;
  cout << endl;

  // test: copy/move assignment
  p_tr = new Kdtree<2,double>(generate_grid(0, 1, 0, -1, 2, 2));
  temp_tr = *p_tr;
  delete p_tr;
  leaves.clear(); temp_tr.report_leaves(leaves);
  cout << "+ copy assignment (1): { ";
  for (const auto &p : leaves) { cout << p << " "; }
  cout << "} (c.f. { (0, -0.5) (0, 0) (0.5, -0.5) (0.5, 0) })" << endl;

  Kdtree<2,double> temp_tr2(generate_grid(0, 1, 0, -1, 2, 2));
  temp_tr = std::move(temp_tr2);
  leaves.clear(); temp_tr.report_leaves(leaves);
  cout << "+ move assignment (1): { ";
  for (const auto &p : leaves) { cout << p << " "; }
  cout << "}, " << temp_tr2.empty();
  cout << " (c.f. { (0, -0.5) (0, 0) (0.5, -0.5) (0.5, 0) }, 1)" << endl;

  cout << endl;

  // test: empty(), size(), get_bounding_box():
  tr1 = generate_grid(0, 1, 0, 1, 100, 100);
  cout << "+ empty() (1): " << tr1.empty() << " (c.f. 0)" << endl;
  cout << "+ size() (1): " << tr1.size() << " (c.f. 10000)" << endl;
  cout << "+ get_bounding_box (1): " << tr1.get_bounding_box(); 
  cout << " (c.f. { (0, 0.99), (0, 0.99) })" << endl;
  cout << endl;

  return 0;
}
