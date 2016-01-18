#include <iostream>
#include <sstream>
#include <vector>
#include <exception>

#include <DefaultNodeAttributes.h>
#include <PointWeights.h>

using namespace std;
using bbrcit::PointWeights;
using bbrcit::DefaultNodeAttributes;

int main() {

  cout << endl;

  using PointAttr = PointWeights<double>;
  using NodeAttr = DefaultNodeAttributes<PointAttr>;

  // test: default constructor
  NodeAttr d0;
  cout << "+ default constructor (1): " << d0 << " (c.f. (1.0) ) " << endl;
  cout << endl;

  // test: one argument constructors
  NodeAttr d1(PointAttr(2.1));
  cout << "+ one argument constructors (1): " << d1 << " (c.f. (2.1) ) " << endl;

  cout << endl;

  // test: set/get_point_attribute
  d0.set_point_attributes(PointAttr(3.2));
  cout << "+ (set|get)_(point|weight) (1): " << d0.get_point_attributes() << " (c.f. (3.2) ) " << endl;
  cout << endl;

  // test: copy/move construct
  d0.set_point_attributes(3.2); 
  NodeAttr d2(d0);
  NodeAttr temp_d0(3.2);
  NodeAttr d3(std::move(temp_d0));
  cout << "+ copy constructor (1): " << d2 << " (c.f. (3.2) ) " << endl;
  cout << "+ move constructor (1): " << d3 << " (c.f. (3.2) ) " << endl;
  cout << endl;

  // test: swap
  d0.set_point_attributes(3.2);
  d1.set_point_attributes(2.3);
  swap(d0, d1);
  cout << "+ swap (1): " << d0 << ", " << d1 << " (c.f. (2.3), (3.2) ) " << endl;

  cout << endl;

  // test: merge
  d0.set_point_attributes(3.2);
  d1.set_point_attributes(2.3);
  d2 = merge(d0, d1); 
  cout << "+ merge (1): " << d2 << " (c.f. (5.5) ) " << endl;

  cout << endl;
  
  return 0;
}
