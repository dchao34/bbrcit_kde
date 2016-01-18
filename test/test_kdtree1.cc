#include <iostream>
#include <vector>
#include <utility>
#include <random>

#include <Rectangle.h>
#include <DecoratedPoint.h>
#include <Kdtree.h>

using namespace std;
using bbrcit::DecoratedPoint;
using bbrcit::Kdtree;

using Kdtree2d = Kdtree<2>;
using Point2d = typename Kdtree2d::PointType;

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
  Kdtree2d tr0;
  cout << "+ Default constructor: " << tr0.empty() << " (cf. 1)" << endl;
  cout << endl;

  // test: vector constructor
  vector<Point2d> points1 = generate_grid(0, 1, 0, 1, 2, 2);
  Kdtree2d tr1(points1);
  vector<Point2d> leaves; tr1.report_leaves(leaves);
  cout << "+ vector constructor (1): " << endl;
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}" << endl;
  cout << "  compare:{ { (0, 0), (1) } { (0, 0.5), (1) } { (0.5, 0), (1) } { (0.5, 0.5), (1) } }" << endl;

  vector<Point2d> points2 = generate_grid(0, -1, 0, -1, 2, 2);
  Kdtree2d tr2(std::move(points2));
  leaves.clear(); tr2.report_leaves(leaves);
  cout << "+ vector constructor (2): " << endl;
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}" << endl;
  cout << "  compare:{ { (-0.5, -0.5), (1) } { (-0.5, 0), (1) } { (0, -0.5), (1) } { (0, 0), (1) } }" << endl;

  cout << endl;

  // test: copy constructor
  Kdtree2d *p_tr = new Kdtree2d(generate_grid(0, 1, 0, -1, 2, 2));
  Kdtree2d tr3(*p_tr);
  delete p_tr;
  leaves.clear(); tr3.report_leaves(leaves);
  cout << "+ copy constructor (1): " << endl;
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}" << endl;
  cout << "  compare:{ { (0, -0.5), (1) } { (0, 0), (1) } { (0.5, -0.5), (1) } { (0.5, 0), (1) } }" << endl;

  Kdtree2d temp_tr(generate_grid(0, -1, 0, 1, 2, 2));
  Kdtree2d tr4(std::move(temp_tr));
  leaves.clear(); tr4.report_leaves(leaves);
  cout << "+ move constructor (1): " << endl;
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}, " << temp_tr.empty() << endl;
  cout << "  compare:{ { (-0.5, 0), (1) } { (-0.5, 0.5), (1) } { (0, 0), (1) } { (0, 0.5), (1) } }, 1" << endl;
  cout << endl;

  // test: swap
  cout << "+ swap (1): " << endl;
  swap(tr3, tr4); 
  leaves.clear(); tr3.report_leaves(leaves);
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}" << endl;
  cout << "  compare:{ { (-0.5, 0), (1) } { (-0.5, 0.5), (1) } { (0, 0), (1) } { (0, 0.5), (1) } }" << endl;

  leaves.clear(); tr4.report_leaves(leaves);
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}" << endl;
  cout << "  compare:{ { (0, -0.5), (1) } { (0, 0), (1) } { (0.5, -0.5), (1) } { (0.5, 0), (1) } }" << endl;

  swap(tr3, tr4); 
  leaves.clear(); tr3.report_leaves(leaves);
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}" << endl;
  cout << "  compare:{ { (0, -0.5), (1) } { (0, 0), (1) } { (0.5, -0.5), (1) } { (0.5, 0), (1) } }" << endl;

  leaves.clear(); tr4.report_leaves(leaves);
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}" << endl;
  cout << "  compare:{ { (-0.5, 0), (1) } { (-0.5, 0.5), (1) } { (0, 0), (1) } { (0, 0.5), (1) } }" << endl;

  cout << endl;

  // test: copy/move assignment
  p_tr = new Kdtree2d(generate_grid(0, 1, 0, -1, 2, 2));
  temp_tr = *p_tr;
  delete p_tr;
  leaves.clear(); temp_tr.report_leaves(leaves);
  cout << "+ copy assignment (1): " << endl;
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}" << endl;
  cout << "  compare:{ { (0, -0.5), (1) } { (0, 0), (1) } { (0.5, -0.5), (1) } { (0.5, 0), (1) } }" << endl;

  Kdtree2d temp_tr2(generate_grid(0, -1, 0, 1, 2, 2));
  temp_tr = std::move(temp_tr2);
  leaves.clear(); temp_tr.report_leaves(leaves);
  cout << "+ move assignment (1): " << endl;
  cout << "  output: { ";for (const auto &p : leaves) { cout << p << " "; } cout << "}, " << temp_tr2.empty() << endl;
  cout << "  compare:{ { (-0.5, 0), (1) } { (-0.5, 0.5), (1) } { (0, 0), (1) } { (0, 0.5), (1) } }, 1" << endl;

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
