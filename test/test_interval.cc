#include <iostream>
#include <utility>
#include <exception>
#include <Interval.h>

using namespace std;
using bbrcit::Interval;

int main() {

  Interval<double> i1(0.1, 0.3);
  Interval<double> i2(0.4, 0.5);
  Interval<double> i3(0.2, 0.39);
  cout << i1 << " " << i2 << " " << i3 << endl;
  cout << i1.length() << " " << i2.length() << " " << i3.length() << endl;
  cout << i1.middle() << " " << i2.middle() << " " << i3.middle() << endl;
  cout << intersect(i1, i2) << " " << intersect(i1, i3) << " " << intersect(i2, i3) << endl;
  swap(i1, i2);
  cout << i1 << " " << i2 << " " << i3 << endl;
  cout << i3.lower() << " " << i3.upper() << endl;
  cout << endl;

  Interval<double> i4({1.1,2.1}), i5({1.2, 2.0});
  cout << i4 << endl;
  cout << i4.contains(1.2) << " " << i4.contains(1.0) << endl;
  cout << i4.contains(i1) << " " << i4.contains(i5) << endl;
  cout << endl;

  Interval<double> i7 = { 3.0, 5.0 };
  cout << i7.min_dist(4.0) << " " << i7.min_dist(1.0) << " " << i7.min_dist(8.0) << endl;
  cout << i7.max_dist(4.0) << " " << i7.max_dist(1.0) << " " << i7.max_dist(8.0) << endl;

  return 0;
}
