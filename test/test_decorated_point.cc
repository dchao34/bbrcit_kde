#include <iostream>
#include <sstream>
#include <vector>
#include <exception>

#include <PointWeights.h>
#include <DecoratedPoint.h>

using namespace std;
using bbrcit::DecoratedPoint;
using bbrcit::PointWeights;
using bbrcit::Point;

int main() {

  using DecoratedPoint2d = DecoratedPoint<2,PointWeights<double>>;

  cout << endl;

  // test: default constructor
  DecoratedPoint2d p0;
  cout << "+ Default constructor: " << p0 << " (c.f. { (0.0, 0.0), (1.0) }) " << endl;
  cout << endl;

  // test: Point<> constructors
  Point<2, double> ep0({1.0, 3.0});
  DecoratedPoint2d p1(ep0);
  cout << "+ DecoratedPoint<> constructors (1): " << p1 << " (c.f. { (1.0, 3.0), (1.0) }) " << endl;

  DecoratedPoint2d p2({2.0, 1.0});
  cout << "+ DecoratedPoint<> constructors (2): " << p2 << " (c.f. { (2.0, 1.0), (1.0) }) " << endl;

  Point<2, double> temp_ep({3.0, 1.0});
  DecoratedPoint2d p3(std::move(temp_ep));
  cout << "+ DecoratedPoint<> constructors (3): " << p3 << " (c.f. { (3.0, 1.0), (1.0) }) " << endl;


  // test: PointT, AttributeT constructors
  ep0 = {1.0, 1.0};
  DecoratedPoint2d p4(ep0, 4.0);
  cout << "+ DecoratedPoint<> constructors (4): " << p4 << " (c.f. { (1.0, 1.0), (4.0) }) " << endl;

  DecoratedPoint2d p5({2.0,2.0}, 3.0);
  cout << "+ DecoratedPoint<> constructors (5): " << p5 << " (c.f. { (2.0, 2.0), (3.0) }) " << endl;

  DecoratedPoint2d::AttributesType attr0(3.0);
  DecoratedPoint2d p6(ep0, attr0);
  cout << "+ DecoratedPoint<> constructors (6): " << p6 << " (c.f. { (1.0, 1.0), (3.0) }) " << endl;

  cout << endl;

  // test: set/get_point, set/get_weight.
  p0.set_point({3.4, 2.2}); p0.set_attributes(0.7);
  cout << "+ (set|get)_(point|attributes) (1): " << p0 << " (c.f. { (3.4, 2.2), (0.7) }) " << endl;
  ep0 = {0, 0}; attr0 = 1.0; p0.set_point(ep0); p0.set_attributes(attr0);
  cout << "+ (set|get)_(point|attributes) (2): " << p0 << " (c.f. { (0.0, 0.0), (1.0) }) " << endl;

  cout << endl;


  // test: copy/move constructors
  p0.set_point({1.0, 2.0}); p0.set_attributes(0.1);
  DecoratedPoint2d temp_p0(p0);
  cout << "+ copy/move constructor (1): " << temp_p0 << " (c.f. { (1.0, 2.0), (0.1) }) " << endl;
  DecoratedPoint2d p7(std::move(temp_p0));
  cout << "+ copy/move constructor (2): " << p7 << " (c.f. { (1.0, 2.0), (0.1) }) " << endl;
  cout << endl;

  // test: copy/move assignment
  p0.set_point({0.0, 0.0}); p0.set_attributes(0.0);
  p1.set_point({0.0, 0.0}); p1.set_attributes(0.0);
  DecoratedPoint2d temp_p1({1.0, 2.0}, 0.1);
  p0 = temp_p1;
  cout << "+ copy/move assignment (1): " << p0 << " (c.f. { (1.0, 2.0), (0.1) }) " << endl;
  p1 = std::move(temp_p1);
  cout << "+ copy/move assignment (2): " << p1 << " (c.f. { (1.0, 2.0), (0.1) }) " << endl;
  cout << endl;


  // test: operator>>
  istringstream sin("0.0 2.0 0.3 abc def 20"); 
  sin >> p0; 
  cout << "+ operator>> (1): " << p0 << " (c.f. { (0.0, 2.0), (0.3) }) " << endl;
  sin >> p0; 
  cout << "+ operator>> (2): " << p0 << " (c.f. { (0.0, 0.0), (1.0) }) " << endl;
  cout << endl;

  // test: operator[]
  p0 = {{ 1.0, 2.0}, 0.3};
  cout << "+ operator[] (1): "; 
  cout << "{ (" << p0[0] << ", " << p0[1] << "), " << p0.attributes() << " } ";
  cout << "(c.f. " << "{ (1.0, 2.0), (0.3) })" << endl;

  const DecoratedPoint2d cp0({1.0, 2.0}, 0.3);
  cout << "+ operator[] (2): "; 
  cout << "{ (" << p0[0] << ", " << p0[1] << "), " << p0.attributes() << " } ";
  cout << "(c.f. " << "{ (1.0, 2.0), (0.3) })" << endl;

  cout << endl;

  // test: swap()
  p0 = {{ 1.0, 2.0}, 0.3}; p1 = {{ 3.0, 2.0}, 0.1};
  swap(p0, p1);
  cout << "+ swap (1): " << endl; 
  cout << "\t " << p0 << " ";
  cout << "(c.f. { (3.0, 2.0), (0.1) })" << endl;
  cout << "\t " << p1 << " ";
  cout << "(c.f. { (1.0, 2.0), (0.3) })" << endl;

  cout << endl;
  
  return 0;
}
