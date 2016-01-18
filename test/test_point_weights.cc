#include <iostream>
#include <sstream>
#include <vector>
#include <exception>

#include <PointWeights.h>

using namespace std;
using bbrcit::PointWeights;

int main() {

  cout << endl;

  using WeightAttDouble = PointWeights<double>;
  using WeightAttInt = PointWeights<int>;

  // test: default constructor
  WeightAttDouble d0;
  cout << "+ default constructor (1): " << d0 << " (c.f. (1.0) ) " << endl;
  WeightAttInt i0;
  cout << "+ default constructor (2): " << d0 << " (c.f. (1) ) " << endl;
  cout << endl;

  // test: one argument constructors
  WeightAttDouble d1(2.1);
  cout << "+ one argument constructors (1): " << d1 << " (c.f. (2.1) ) " << endl;

  WeightAttInt i1(2);
  cout << "+ one argument constructors (2): " << i1 << " (c.f. (2) ) " << endl;

  cout << endl;

  // test: set/get_weight
  d0.set_weight(3.2); i0.set_weight(3);
  cout << "+ (set|get)_(point|weight) (1): " << d0 << " (c.f. (3.2) ) " << endl;
  cout << "+ (set|get)_(point|weight) (2): " << i0 << " (c.f. (3) ) " << endl;
  cout << endl;

  // test: copy/move construct
  d0.set_weight(3.2); i0.set_weight(3);
  WeightAttDouble d2(d0);
  WeightAttInt i2(i0);
  WeightAttDouble temp_d0(3.2);
  WeightAttInt temp_i0(3);
  WeightAttDouble d3(std::move(temp_d0));
  WeightAttInt i3(std::move(temp_i0));
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

  // test: merge
  d0.set_weight(3.2); i0.set_weight(3);
  d1.set_weight(2.3); i1.set_weight(2);
  d2 = merge(d0, d1); i2 = merge(i0, i1);
  cout << "+ merge (1): " << d2 << " (c.f. (5.5) ) " << endl;
  cout << "+ merge (2): " << i2 << " (c.f. (5) ) " << endl;

  cout << endl;
  
  return 0;
}
