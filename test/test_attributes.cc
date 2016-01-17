#include <iostream>
#include <sstream>
#include <vector>
#include <exception>

#include <PointAttributes.h>

using namespace std;
using bbrcit::MinimalAttributes;

int main() {

  cout << endl;

  using MinAttDouble = MinimalAttributes<double>;
  using MinAttInt = MinimalAttributes<int>;

  // test: default constructor
  MinAttDouble d0;
  cout << "+ default constructor (1): " << d0 << " (c.f. (1.0) ) " << endl;
  MinAttInt i0;
  cout << "+ default constructor (2): " << d0 << " (c.f. (1) ) " << endl;
  cout << endl;

  // test: one argument constructors
  MinAttDouble d1(2.1);
  cout << "+ one argument constructors (1): " << d1 << " (c.f. (2.1) ) " << endl;

  MinAttInt i1(2);
  cout << "+ one argument constructors (2): " << i1 << " (c.f. (2) ) " << endl;

  cout << endl;

  // test: set/get_weight
  d0.set_weight(3.2); i0.set_weight(3);
  cout << "+ (set|get)_(point|weight) (1): " << d0 << " (c.f. (3.2) ) " << endl;
  cout << "+ (set|get)_(point|weight) (2): " << i0 << " (c.f. (3) ) " << endl;
  cout << endl;

  // test: copy/move construct
  d0.set_weight(3.2); i0.set_weight(3);
  MinAttDouble d2(d0);
  MinAttInt i2(i0);
  MinAttDouble temp_d0(3.2);
  MinAttInt temp_i0(3);
  MinAttDouble d3(std::move(temp_d0));
  MinAttInt i3(std::move(temp_i0));
  cout << "+ copy constructor (1): " << d2 << " (c.f. (3.2) ) " << endl;
  cout << "+ copy constructor (2): " << i2 << " (c.f. (3) ) " << endl;
  cout << "+ move constructor (1): " << d3 << " (c.f. (3.2) ) " << endl;
  cout << "+ move constructor (2): " << i3 << " (c.f. (3) ) " << endl;
  cout << endl;

  // test: swap
  d0.set_weight(3.2); i0.set_weight(3);
  d1.set_weight(2.3); i1.set_weight(2);
  swap(d0, d1); swap(i0, i1);
  cout << "+ copy constructor (1): " << d0 << ", " << d1 << " (c.f. (2.3), (3.2) ) " << endl;
  cout << "+ copy constructor (2): " << i0 << ", " << i1 << " (c.f. (2), (3) ) " << endl;

  cout << endl;

  // test: combine_weights
  d0.set_weight(3.2); i0.set_weight(3);
  d1.set_weight(2.3); i1.set_weight(2);
  d2 = add_weights(d0, d1); i2 = add_weights(i0, i1);
  cout << "+ add_weights (1): " << d2 << " (c.f. (5.5) ) " << endl;
  cout << "+ add_weights (2): " << i2 << " (c.f. (5) ) " << endl;

  cout << endl;
  
  return 0;
}
