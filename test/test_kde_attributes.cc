#include <iostream>
#include <sstream>
#include <vector>
#include <exception>

#include <KdeAttributes.h>

using namespace std;
using bbrcit::KdeAttributes;

int main() {

  cout << endl;

  using PointAttr = KdeAttributes<double>;

  // test: default constructor
  PointAttr d0;
  cout << "+ default constructor (1): " << d0 << " (c.f. (1.0, 0, 0) " << endl;
  cout << endl;

  // test: one argument constructors
  PointAttr d1(2.1, 1.0, 10.0);
  cout << "+ one argument constructors (1): " << d1 << " (c.f. (2.1, 1.0, 10.0) ) " << endl;

  cout << endl;

  // test: set/get attributes
  d0.set_weight(3.2); d0.set_lower(1.0); d0.set_upper(10.0);
  cout << "+ set/get attributes (1): " << d0.weight() << " (c.f. (3.2) ) " << endl;
  cout << "+ set/get attributes (1): " << d0.lower() << " (c.f. (1.0) ) " << endl;
  cout << "+ set/get attributes (1): " << d0.upper() << " (c.f. (10.0) ) " << endl;
  cout << "+ set/get attributes (1): " << d0.middle() << " (c.f. (5.5) ) " << endl;
  cout << endl;

  // test: copy/move construct
  d0.set_weight(3.2); d0.set_lower(1.0); d0.set_upper(10.0);
  PointAttr d2(d0);
  PointAttr temp_d0(3.2, 1.0, 10.0);
  PointAttr d3(std::move(temp_d0));
  cout << "+ copy constructor (1): " << d2 << " (c.f. (3.2, 1.0, 10.0) ) " << endl;
  cout << "+ move constructor (1): " << d3 << " (c.f. (3.2, 1.0, 10.0) ) " << endl;
  cout << endl;

  // test: swap
  d0.set_weight(3.2); d0.set_lower(1.0); d0.set_upper(10.0);
  d1.set_weight(2.3); d1.set_lower(2.0); d1.set_upper(9.0);
  swap(d0, d1);
  cout << "+ swap (1): " << d0 << ", " << d1 << " (c.f. (2.3, 2.0, 9.0), (3.2, 1.0, 10.0) ) " << endl;

  cout << endl;

  // test: merge
  d0.set_weight(3.2); d0.set_lower(1.0); d0.set_upper(10.0);
  d1.set_weight(2.3); d1.set_lower(2.0); d1.set_upper(9.0);
  d2 = merge(d0, d1); 
  cout << "+ merge (1): " << d2 << " (c.f. (5.5, 1.0, 10.0) ) " << endl;

  cout << endl;
  
  return 0;
}
