#include <iostream>
#include <utility>
#include <Rectangle.h>
#include <Point.h>
#include <DecoratedPoint.h>

using namespace std;
using bbrcit::Point;
using bbrcit::Rectangle;
using bbrcit::DecoratedPoint;

int main() {

  cout << endl;

  // test: default constructor
  Rectangle<2> r0;
  cout << "+ Default constructor: " << r0 << " (cf. { (0.0, 0.0), (0.0, 0.0) }) " << endl;
  cout << endl;

  // test: two argument constructor
  Rectangle<2> r1({0.0, 0.0}, {1.0, 1.0});
  cout << "+ Two argument constructor: " << r1 << " (cf. { (0.0, 1.0), (0.0, 1.0) }) " << endl;
  cout << endl;

  // test: Copy/move constructor
  Rectangle<2> r2(r1);
  cout << "+ Copy/move constructor (1): " << r2 << " (cf. { (0.0, 1.0), (0.0, 1.0) }) " << endl;
  Rectangle<2> temp_r1({0.0,0.0}, {1.0,1.0});
  Rectangle<2> r3(std::move(temp_r1));
  cout << "+ Copy/move constructor (2): " << r3 << " (cf. { (0.0, 1.0), (0.0, 1.0) }) " << endl;
  cout << endl;

  // test: Copy/move assignment
  Rectangle<2> r4, r5({0.0, 0.0}, {1.0,1.0});
  r4 = r5;
  cout << "+ Copy/move assignment (1): " << r4 << " (cf. { (0.0, 1.0), (0.0, 1.0) }) " << endl;
  Rectangle<2> temp_r2({0.0, 0.0}, {2.0,2.0});
  r5 = std::move(temp_r2);
  cout << "+ Copy/move assignment (2): " << r5 << " (cf. { (0.0, 2.0), (0.0, 2.0) }) " << endl;
  cout << endl;

  // test: operator[]
  const Rectangle<2> cp0({2.0, 0.0}, {4.0, 3.0});
  cout << "+ operator[] (1): " << "{ "<< r1[0] << ", " << r1[1] << " }" << " (cf. { (0.0, 1.0), (0.0, 1.0) }) " << endl;
  cout << "+ operator[] (2): " << "{ "<< cp0[0] << ", " << cp0[1] << " }" << " (cf. { (2.0, 4.0), (0.0, 3.0) }) " << endl;
  cout << endl;

  // test: resize()
  Rectangle<2> r6({0,0}, {1,1});
  r6.resize(0, {0, 2}); r6.resize(1, {0, 2});
  cout << "+ resize() (1): " << r6 << " (cf. { (0.0, 2.0), (0.0, 2.0) }) " << endl;
  r6.resize({1,1}, {3,3});
  cout << "+ resize() (2): " << r6 << " (cf. { (1.0, 3.0), (1.0, 3.0) }) " << endl;
  cout << endl;

  // test: lower_halfspace(), upper_halfspace()
  r0.resize({0,0}, {4,4});
  r1 = r0.lower_halfspace(0, 1); r2 = r0.upper_halfspace(0, 1);
  cout << "+ lower_halfspace() (1): " << r1 << " (cf. { (0.0, 1.0), (0.0, 4.0) }) " << endl;
  cout << "+ upper_halfspace() (1): " << r2 << " (cf. { (1.0, 4.0), (0.0, 4.0) }) " << endl;
  r1 = r0.lower_halfspace(1, 3); r2 = r0.upper_halfspace(1, 3);
  cout << "+ lower_halfspace() (1): " << r1 << " (cf. { (0.0, 4.0), (0.0, 3.0) }) " << endl;
  cout << "+ upper_halfspace() (1): " << r2 << " (cf. { (0.0, 4.0), (3.0, 4.0) }) " << endl;
  cout << endl;

  // test: swap()
  r0.resize({1,1}, {3,3}); r1.resize({0,0}, {1,1});
  swap(r0, r1);
  cout << "+ swap() (1): " << r0 << ", " << r1;
  cout << " (cf. { (0.0, 1.0), (0.0, 1.0) }, { (1.0, 3.0), (1.0, 3.0) }) " << endl;
  cout << endl;

  // test: contains()
  Point<2> p0, p1({-3,0}), p2({0, 3});
  r0.resize({-2,-2}, {2,2}); 
  cout << "+ contains() (1): "; 
  cout << r0.contains(p0) << ", ";
  cout << r0.contains(p1) << ", ";
  cout << r0.contains(p2);
  cout << " (cf. 1, 0, 0)" << endl;

  r1.resize({-1,-1}, {1,1}); r2.resize({-2,-1}, {2,1}); r3.resize({-3,-1}, {2,1}); r4.resize({3,0}, {4,1});
  cout << "+ contains() (2): "; 
  cout << r0.contains(r0) << ", ";
  cout << r0.contains(r1) << ", ";
  cout << r0.contains(r2) << ", ";
  cout << r0.contains(r3) << ", ";
  cout << r0.contains(r4);
  cout << " (cf. 1, 1, 1, 0, 0)" << endl;

  DecoratedPoint<2> dp0, dp1({-3,0}), dp2({0, 3});
  r0.resize({-2,-2}, {2,2}); 
  cout << "+ contains() (3): "; 
  cout << r0.contains(dp0) << ", ";
  cout << r0.contains(dp1) << ", ";
  cout << r0.contains(dp2);
  cout << " (cf. 1, 0, 0)" << endl;
  cout << endl;

  // test: intersect()
  r0.resize({-2,-2}, {2,2}); 
  r1.resize({-1,-1}, {1,1}); r2.resize({-2,-3}, {2,-1}); r3.resize({3,-2}, {5,3});
  cout << "+ intersect() (1): "; 
  cout << intersect(r0, r0) << ", ";
  cout << intersect(r0, r1) << ", ";
  cout << intersect(r0, r2) << ", ";
  cout << intersect(r0, r3);
  cout << " (cf. 1, 1, 1, 0)" << endl;
  cout << endl;

  // test: min_dist(), max_dist()
  p0.reset({0,0});
  r1.resize({2,0}, {4,3}); 
  cout << "+ min_dist(), max_dist() (1): " << r1.min_dist(p0) << ", " << r1.max_dist(p0) << " (c.f 2, 5)" << endl; 
  r1.resize({1,-3}, {4,-1}); 
  cout << "+ min_dist(), max_dist() (2): " << r1.min_dist(p0) << ", " << r1.max_dist(p0) << " (c.f 1.41421, 5)" << endl; 
  r1.resize({-1,-1}, {1,1}); 
  cout << "+ min_dist(), max_dist() (3): " << r1.min_dist(p0) << ", " << r1.max_dist(p0) << " (c.f 0, 1.41421)" << endl; 
  r1.resize({0,0}, {1,1}); 
  cout << "+ min_dist(), max_dist() (4): " << r1.min_dist(p0) << ", " << r1.max_dist(p0) << " (c.f 0, 1.41421)" << endl; 
  r1.resize({-0.5,-0.5}, {1,1}); 
  cout << "+ min_dist(), max_dist() (5): " << r1.min_dist(p0) << ", " << r1.max_dist(p0) << " (c.f 0, 1.41421)" << endl; 
  cout << endl;


  return 0;
}
