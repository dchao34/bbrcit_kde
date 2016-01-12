#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <exception>
#include <utility>
#include <Point.h>

using namespace std;
using bbrcit::Point;

Point<3> f() {
  return Point<3>{1,2,3};
}

int main() {

  cout << endl;

  Point<2> p0, p1, p2, p3, p4, p5, p6;
  const Point<2> cp0({0.0,0.0}), cp1({1.0,1.0}), cp2({2.0,1.0});

  // test: default constructor
  cout << "+ Default constructor gives: " << p0 << " (should be 0.0, 0.0). " << endl;
  cout << endl;

  // test: initializer list constructor
  Point<2> test_p1({0.0, 1.0});
  Point<2> test_p2 = {0.0, 2.0};
  Point<2> test_p8({1.0});
  cout << "+ Initilizer list constructor (1): " << test_p1 << " (should be (0.0, 1.0)). " << endl;
  cout << "+ Initilizer list constructor (2): " << test_p2 << " (should be (0.0, 2.0)). " << endl;
  cout << "+ Initilizer list constructor (3): " << test_p8 << " (should be (1.0, 0.0)). " << endl;
  cout << endl;

  // test: copy/move constructor
  Point<2> test_p3 = { 0.0, 1.0 };
  Point<2> test_p4 = { 0.0, 2.0 };
  Point<2> test_p5(test_p3);
  Point<2> test_p6(std::move(test_p4));
  cout << "+ Copy constructor: " << test_p5 << ", " << test_p3 << " (should both be (0.0, 1.0)). " << endl;
  cout << "+ Move constructor: " << test_p6 << " (should be (0.0, 2.0)). " << endl;
  cout << endl;

  // test: copy/move assignment
  p0 = { 0.0, 1.0 };
  cout << "+ Copy/Move assignment (1): " << p0 << " (should be (0.0, 1.0)). " << endl;
  p0 = p0;
  cout << "+ Copy/Move assignment (2): " << p0 << " (should be (0.0, 1.0)). " << endl;
  p1 = p0;
  cout << "+ Copy/Move assignment (3): " << p1 << " (should be (0.0, 1.0)). " << endl;
  Point<2> test_p7 = { 1.0, 2.0 };
  p0 = std::move(test_p7);
  cout << "+ Copy/Move assignment (4): " << p0 << " (should be (1.0, 2.0)). " << endl;
  cout << endl;

  // test: swap()
  p0 = { 0.0, 1.0 }; p1 = { 1.0, 0.0 };
  swap(p0, p0);
  cout << "+ swap(): " << p0 << " " << p1 << " (should be (0.0, 1.0), (1.0, 0.0). " << endl;
  swap(p0, p1);
  cout << "+ swap(): " << p0 << " " << p1 << " (should be (1.0, 0.0), (0.0, 1.0). " << endl;
  cout << endl;

  // test: reset().
  p0.reset({1.0, 1.0});
  cout << "+ reset() (1): " << p0 << " (should be (1.0, 1.0)). " << endl;
  p0.reset({1.0});
  cout << "+ reset() (2): " << p0 << " (should be (1.0, 0.0)). " << endl;
  p0.reset({});
  cout << "+ reset() (3): " << p0 << " (should be (0.0, 0.0)). " << endl;
  cout << endl; 

  // test: operator>>
  p0.reset({});
  istringstream sin("3.0 2.0 1.0 abc 3.0 3.0");
  sin >> p0;
  cout << "+ operator>>() (1): " << p0 << " (should be (3.0, 2.0)). " << endl;
  sin >> p0;
  cout << "+ operator>>() (2): " << p0 << " (should be (0.0, 0.0)). " << endl;
  sin.clear();
  sin >> p0;
  cout << "+ operator>>() (3): " << p0 << " (should be (3.0, 3.0)). " << endl;
  sin >> p0;
  cout << "+ operator>>() (3): " << p0 << " (should be (0.0, 0.0)). " << endl;
  sin.clear();
  cout << endl; 

  // test: operator[]
  p0.reset({1.0, 2.0});
  cout << "+ operator[] (1): " << p0[0] << ", " << p0[1] << " (should be 1.0, 2.0). " << endl;
  cout << "+ operator[] (2): " << cp2[0] << ", " << cp2[1] << " (should be 2.0, 1.0). " << endl;
  p0[0] = 3.0; p0[1] = 3.0; 
  cout << "+ operator[] (3): " << p0[0] << ", " << p0[1] << " (should be 3.0, 3.0). " << endl;
  cout << endl; 

  // test: operator*=
  p0.reset({1.0, 2.0});
  p0 *= 3.0;
  cout << "+ operator*= (1): " << p0 << " (should be (3.0, 6.0)). " << endl;
  cout << endl; 

  // test: operator/=
  p0.reset({3.0, 6.0});
  p0 /= 3.0;
  cout << "+ operator/= (1): " << p0 << " (should be (1.0, 2.0)). " << endl;
  try { 
    p0 /= 0.0;
    cout << "+ operator/= (2): ERROR: mishandled division by zero. " << endl;
  } catch (exception &e) {
    cout << "+ operator/= (2): passed: division by zero. " << endl;
  }
  cout << endl; 

  // test: operator+=
  p0.reset({1.0, 2.0}); p1.reset({2.0, 1.0});
  p0 += p1;
  cout << "+ operator+= (1): " << p0 << " (should be (3.0, 3.0)). " << endl;
  p0 += p0;
  cout << "+ operator+= (2): " << p0 << " (should be (6.0, 6.0)). " << endl;
  cout << endl; 

  // test: operator-=
  p0.reset({3.0, 2.0}); p1.reset({1.0, 1.0});
  p0 -= p1;
  cout << "+ operator-= (1): " << p0 << " (should be (2.0, 1.0)). " << endl;
  p0 -= p0;
  cout << "+ operator-= (2): " << p0 << " (should be (0.0, 0.0)). " << endl;
  cout << endl; 

  // test: operator*
  p0.reset({1.0, 2.0});
  p1 = p0 * 2;
  cout << "+ operator* (1): " << p1 << " (should be (2.0, 4.0)). " << endl;
  p1 = p1 * 2;
  cout << "+ operator* (2): " << p1 << " (should be (4.0, 8.0)). " << endl;
  p1 = 0.5 * p1;
  cout << "+ operator* (3): " << p1 << " (should be (2.0, 4.0)). " << endl;
  cout << endl; 

  // test: operator/
  p0.reset({4.0, 8.0});
  p1 = p0 / 2.0;
  cout << "+ operator/ (1): " << p1 << " (should be (2.0, 4.0)). " << endl;
  p1 = p1 / 2;
  cout << "+ operator/ (2): " << p1 << " (should be (1.0, 2.0)). " << endl;
  cout << endl; 

  // test: operator+
  p0.reset({1.0, 1.0}); p1.reset({-1.0, -1.0});
  p2 = p0 + p0;
  cout << "+ operator+ (1): " << p2 << " (should be (2.0, 2.0)). " << endl;
  p2 = p2 + p1;
  cout << "+ operator+ (1): " << p2 << " (should be (1.0, 1.0)). " << endl;
  cout << endl; 

  // test: operator-
  p0.reset({1.0, 1.0}); p1.reset({-1.0, -1.0});
  p2 = p1 - p0;
  cout << "+ operator- (1): " << p2 << " (should be (-2.0, -2.0)). " << endl;
  p2 = p2 - p2;
  cout << "+ operator- (2): " << p2 << " (should be (0.0, 0.0)). " << endl;
  cout << endl; 


  return 0;
}
