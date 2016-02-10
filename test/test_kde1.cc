#include <iostream>
#include <vector>
#include <utility>
#include <random>

#include <Rectangle.h>
#include <DecoratedPoint.h>
#include <KernelDensity.h>

#include "kde_test_utils.h"

using namespace std;
using bbrcit::DecoratedPoint;
using bbrcit::KernelDensity;

using Kde2d = KernelDensity<2>;
using DataPoint2d = typename Kde2d::DataPointType;
using KernelType = typename KernelDensity<2>::KernelType;

int main() {

  cout << endl;

  // test: default constructor
  Kde2d tr0;
  cout << "+ Default constructor: " << tr0.kernel().bandwidth() << ", " << tr0.size() << " (cf. 1, 0)" << endl;
  cout << endl;

  // test: vector constructor
  vector<DataPoint2d> points; 
  generate_2dgrid(points, 0, 1, 2, 0, 1, 2);
  Kde2d tr1(points); tr1.kernel().set_bandwidth(1);
  cout << "+ vector constructor (1): " << endl;
  cout << "  output: " << endl;
  for (const auto &p : tr1.points()) { cout << "\t" << p << endl; }
  cout << "\t" << tr1.kernel().bandwidth() << endl;
  cout << "  compare: " << endl;
  cout << "\t" << "{ (0, 0), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0, 0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0.5, 0), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0.5, 0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "1" << endl;
  cout << endl;

  generate_2dgrid(points, -1, 0, 2, -1, 0, 2);
  Kde2d tr2(std::move(points)); tr2.kernel().set_bandwidth(2);
  cout << "+ vector constructor (2): " << endl;
  cout << "  output: " << endl;
  for (const auto &p : tr2.points()) { cout << "\t" << p << endl; }
  cout << "\t" << tr2.kernel().bandwidth() << endl;
  cout << "  compare: " << endl;
  cout << "\t" << "{ (-1, -1), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-1, -0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-0.5, -1), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-0.5, -0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "2" << endl;
  cout << endl;

  cout << endl;

  // test: copy constructor
  generate_2dgrid(points, 0, 1, 2, 0, 1, 2);
  Kde2d *p_tr = new Kde2d(points);
  Kde2d tr3(*p_tr);
  delete p_tr;
  cout << "+ copy constructor (1): " << endl;
  cout << "  output: " << endl;
  for (const auto &p : tr3.points()) { cout << "\t" << p << endl; }
  cout << "\t" << tr3.kernel().bandwidth() << endl;
  cout << "  compare: " << endl;
  cout << "\t" << "{ (0, 0), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0, 0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0.5, 0), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0.5, 0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "1" << endl;
  cout << endl;

  generate_2dgrid(points, -1, 0, 2, -1, 0, 2);
  Kde2d temp_tr(points); temp_tr.kernel().set_bandwidth(2);
  Kde2d tr4(std::move(temp_tr));
  cout << "+ move constructor (1): " << endl;
  cout << "  output: " << endl;
  for (const auto &p : tr4.points()) { cout << "\t" << p << endl; }
  cout << "\t" << tr4.kernel().bandwidth() << endl;
  cout << "  compare: " << endl;
  cout << "\t" << "{ (-1, -1), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-1, -0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-0.5, -1), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-0.5, -0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "2" << endl;
  cout << endl;

  // test: swap
  cout << "+ swap (1): " << endl;
  cout << "  output: " << endl;
  swap(tr3, tr4); 
  for (const auto &p : tr3.points()) { cout << "\t" << p << endl; }
  cout << "\t" << tr3.kernel().bandwidth() << endl;
  for (const auto &p : tr4.points()) { cout << "\t" << p << endl; }
  cout << "\t" << tr4.kernel().bandwidth() << endl;
  cout << "  compare: " << endl;
  cout << "\t" << "{ (-1, -1), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-1, -0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-0.5, -1), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-0.5, -0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "2" << endl;
  cout << "\t" << "{ (0, 0), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0, 0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0.5, 0), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0.5, 0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "1" << endl;

  cout << endl;

  // test: copy/move assignment
  generate_2dgrid(points, 0, 1, 2, 0, 1, 2);
  p_tr = new Kde2d(points, 2);
  temp_tr = *p_tr;
  delete p_tr;
  cout << "+ copy assignment (1): " << endl;
  cout << "  output: " << endl;
  for (const auto &p : temp_tr.points()) { cout << "\t" << p << endl; }
  cout << "\t" << temp_tr.kernel().bandwidth() << endl;
  cout << "  compare: " << endl;
  cout << "\t" << "{ (0, 0), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0, 0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0.5, 0), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (0.5, 0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "1" << endl;
  cout << endl;

  generate_2dgrid(points, -1, 0, 2, -1, 0, 2);
  Kde2d temp_tr2(points, 2);
  temp_tr = std::move(temp_tr2);
  cout << "+ move assignment (1): " << endl;
  for (const auto &p : temp_tr.points()) { cout << "\t" << p << endl; }
  cout << "\t" << tr4.kernel().bandwidth() << endl;
  cout << "  compare: " << endl;
  cout << "\t" << "{ (-1, -1), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-1, -0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-0.5, -1), (0.25, 0, 0) }" << endl;
  cout << "\t" << "{ (-0.5, -0.5), (0.25, 0, 0) }" << endl;
  cout << "\t" << "2" << endl;
  cout << endl;

  cout << endl;

  // test: size(), bandwidth():
  generate_2dgrid(points, 0, 1, 100, 0, 1, 100);
  tr1 = Kde2d(points);
  tr1.set_kernel(KernelType(3));
  cout << "+ size() (1): " << tr1.size() << " (c.f. 10000)" << endl;
  cout << "+ bandwidth() (1): " << tr1.kernel().bandwidth() << " (c.f. 3)" << endl;
  cout << endl;


  return 0;
}
