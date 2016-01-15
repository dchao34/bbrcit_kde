#include <iostream>
#include <sstream>
#include <vector>
#include <exception>

#include <PointW.h>

using namespace std;
using bbrcit::PointW;
using bbrcit::Point;

int main() {

  cout << endl;

  // test: default constructor
  PointW<2,double> p0;
  cout << "+ Default constructor: " << p0 << " (c.f. { (0.0, 0.0), 1.0 }) " << endl;
  cout << endl;

  // test: Point<> constructors
  Point<2, double> ep0({1.0, 3.0});
  PointW<2,double> p1(ep0);
  cout << "+ PointW<> constructors (1): " << p1 << " (c.f. { (1.0, 3.0), 1.0 }) " << endl;

  PointW<2,double> p2({2.0, 1.0});
  cout << "+ PointW<> constructors (2): " << p2 << " (c.f. { (2.0, 1.0), 1.0 }) " << endl;

  Point<2, double> temp_ep({3.0, 1.0});
  PointW<2,double> p3(std::move(temp_ep));
  cout << "+ PointW<> constructors (3): " << p3 << " (c.f. { (3.0, 1.0), 1.0 }) " << endl;

  // test: Point<>, WeightT constructors
  ep0 = {1.0, 1.0};
  PointW<2,double> p4(ep0, 4.0);
  cout << "+ PointW<> constructors (4): " << p4 << " (c.f. { (1.0, 1.0), 4.0 }) " << endl;

  PointW<2,double> p5({2.0,2.0}, 3.0);
  cout << "+ PointW<> constructors (5): " << p5 << " (c.f. { (2.0, 2.0), 3.0 }) " << endl;

  cout << endl;

  // test: set/get_point, set/get_weight.
  p0.set_point({3.4, 2.2}); p0.set_weight(0.7);
  cout << "+ (set|get)_(point|weight) (1): " << p0 << " (c.f. { (3.4, 2.2), 0.7 }) " << endl;
  ep0 = {0, 0}; p0.set_point(ep0); p0.set_weight(1.0);
  cout << "+ (set|get)_(point|weight) (2): " << p0 << " (c.f. { (0.0, 0.0), 1.0 }) " << endl;

  cout << endl;

  // test: copy/move constructors
  p0.set_point({1.0, 2.0}); p0.set_weight(0.1);
  PointW<2,double> temp_p0(p0);
  cout << "+ copy/move constructor (1): " << temp_p0 << " (c.f. { (1.0, 2.0), 0.1 }) " << endl;
  PointW<2,double> p6(std::move(temp_p0));
  cout << "+ copy/move constructor (2): " << p6 << " (c.f. { (1.0, 2.0), 0.1 }) " << endl;
  cout << endl;

  // test: copy/move assignment
  p0.set_point({0.0, 0.0}); p0.set_weight(0.0);
  p1.set_point({0.0, 0.0}); p1.set_weight(0.0);
  PointW<2,double> temp_p1({1.0, 2.0}, 0.1);
  p0 = temp_p1;
  cout << "+ copy/move assignment (1): " << p0 << " (c.f. { (1.0, 2.0), 0.1 }) " << endl;
  p1 = std::move(temp_p1);
  cout << "+ copy/move assignment (2): " << p1 << " (c.f. { (1.0, 2.0), 0.1 }) " << endl;
  cout << endl;

  // test: operator>>
  istringstream sin("0.0 2.0 0.3 abc def 20"); 
  sin >> p0; 
  cout << "+ operator>> (1): " << p0 << " (c.f. { (0.0, 2.0), 0.3 }) " << endl;
  sin >> p0; 
  cout << "+ operator>> (2): " << p0 << " (c.f. { (0.0, 0.0), 1.0 }) " << endl;
  cout << endl;

  // test: operator[]
  p0 = {{ 1.0, 2.0}, 0.3};
  cout << "+ operator[] (1): "; 
  cout << "{ (" << p0[0] << ", " << p0[1] << "), " << p0.get_weight() << " } ";
  cout << "(c.f. " << "{ (1.0, 2.0), 0.3 })" << endl;

  const PointW<2,double> cp0({1.0, 2.0}, 0.3);
  cout << "+ operator[] (2): "; 
  cout << "{ (" << p0[0] << ", " << p0[1] << "), " << p0.get_weight() << " } ";
  cout << "(c.f. " << "{ (1.0, 2.0), 0.3 })" << endl;

  cout << endl;
  
  return 0;
}
