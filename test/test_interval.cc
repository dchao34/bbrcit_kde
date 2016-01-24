#include <iostream>
#include <utility>
#include <exception>
#include <Interval.h>

using namespace std;
using bbrcit::Interval;

int main() {

  cout << endl;

  Interval<double> i0, i1, i2, i3, i4, i5, i6; 

  // test: default constructor
  cout << "+ Default constructor gives: " << i0 << " (should be (0.0, 0.0)). " << endl;
  cout << endl;

  // test: two argument constructor (1)
  Interval<double> test_i1(0.0, 1.0);
  cout << "+ Two argument constructor: " << test_i1 << " (should be (0.0, 1.0)). " << endl;
  cout << endl;

  // test: two argument constructor (2)
  try {
    Interval<double> bad1(2.0, 0.0);
    cout << "+ Two argument constructor: ERROR: violates invariant. " << endl;
  } catch (exception &e) {
    cout << "+ Two argument constructor: respects invariant. " << endl;
  }
  cout << endl;

  // test: resize(const T&, const T&) (1)
  i3.resize(0.0, 2.0);
  cout << "+ resize(): " << i3 << " (should be (0.0, 2.0)). " << endl;
  cout << endl;

  // test: resize(const T&, const T&) (2)
  try {
    i3.resize(1.0, 0.0);
    cout << "+ resize(): ERROR: violates invariant. " << endl;
  } catch (exception &e) {
    cout << "+ resize(): " << i3 << " (should still be (0.0, 2.0)). " << endl;
  }
  cout << endl;

  // test: lower(), middle(), upper() (1)
  i1.resize(1.0, 3.0);
  cout << "+ lower(), middle(), upper(), length(): ";
  cout << i1.lower() << ", " << i1.middle() << ", " << i1.upper() << ", " << i1.length();
  cout << " (should be: 1.0, 2.0, 3.0, 2.0)" << endl;
  cout << endl;

  // test: contains(const T&)
  i1.resize(0.0, 1.0); 
  cout << "+ contains(const T&): ";
  cout << i1.contains(-1.0) << ", " << i1.contains(0.5) << ", " << i1.contains(1.5);
  cout << " (should be 0, 1, 0)" << endl;
  cout << endl;

  // test: contains(const Interval<T>&):
  i1.resize(-1.0, 1.0); i2.resize(-0.5, 0.5); i3.resize(1.5, 2.0); i4.resize(0.0, 1.1); i5.resize(0.0, 1.0);
  cout << "+ contains(const Interval<T>&): ";
  cout << i1.contains(i1) << ", ";
  cout << i1.contains(i2) << ", ";
  cout << i1.contains(i3) << ", ";
  cout << i1.contains(i4) << ", ";
  cout << i1.contains(i5);
  cout << " (should be 1, 1, 0, 0, 1)" << endl;
  cout << endl;

  // test: intersect(const Interval<T>&, const Interval<T>&):
  cout << "+ intersect(const Interval<T>&, const Interval<T>&): ";
  i1.resize(-1.0, 1.0); i2.resize(-1.5, -1.1); i3.resize(-1.5, -0.5); 
  i4.resize(-0.5, 0.5); i5.resize(0.5, 1.1); i6.resize(1.1, 1.5);
  cout << intersect(i1, i1) << ", ";
  cout << intersect(i1, i2) << ", ";
  cout << intersect(i1, i3) << ", ";
  cout << intersect(i1, i4) << ", ";
  cout << intersect(i1, i5) << ", ";
  cout << intersect(i1, i6);
  cout << " (should be 1, 0, 1, 1, 1, 0)" << endl;
  cout << endl;

  // test: min_dist(const &T): 
  cout << "+ min_dist(const &T): ";
  i1.resize(1.0, 3.0);
  cout << i1.min_dist(-1) << ", ";
  cout << i1.min_dist(1.5) << ", ";
  cout << i1.min_dist(3.5);
  cout << " (should be 2, 0, 0.5)" << endl;
  cout << endl;

  // test: max_dist(const &T): 
  cout << "+ max_dist(const &T): ";
  i1.resize(1.0, 4.0);
  cout << i1.max_dist(-1) << ", ";
  cout << i1.max_dist(2) << ", ";
  cout << i1.max_dist(3) << ", ";
  cout << i1.max_dist(5);
  cout << " (should be 5, 2, 2, 4)" << endl;
  cout << endl;

  // test: min_dist(const Interval<T>&): 
  cout << "+ min_dist(const Interval<T>&): ";
  i1.resize(10.0, 20.0);
  i2.resize(8, 9); cout << i1.min_dist(i2) << ", ";
  i2.resize(8, 12); cout << i1.min_dist(i2) << ", ";
  i2.resize(14, 16); cout << i1.min_dist(i2) << ", ";
  i2.resize(17, 23); cout << i1.min_dist(i2) << ", ";
  i2.resize(23, 30); cout << i1.min_dist(i2);
  cout << " (c.f. 1, 0, 0, 0, 3)" << endl;
  cout << endl;

  // test: max_dist(const Interval<T>&): 
  cout << "+ max_dist(const Interval<T>&): ";
  i1.resize(10.0, 20.0);
  i2.resize(8, 9); cout << i1.max_dist(i2) << ", ";
  i2.resize(8, 12); cout << i1.max_dist(i2) << ", ";
  i2.resize(14, 17); cout << i1.max_dist(i2) << ", ";
  i2.resize(17, 24); cout << i1.max_dist(i2) << ", ";
  i2.resize(23, 30); cout << i1.max_dist(i2);
  cout << " (c.f. 12, 12, 7, 14, 20)" << endl;
  cout << endl;

  return 0;
}
