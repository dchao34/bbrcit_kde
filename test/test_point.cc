#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <exception>
#include <utility>
#include <Point.h>

using namespace std;

Point<3> f() {
  return Point<3>{1,2,3};
}

int main() {
  Point<3> p1, p2 = { 1, 2.0, 3.2 };
  cout << p1 << ", "; cout << p2 << endl;
  swap(p1, p2);
  cout << p1 << ", "; cout << p2 << endl;
  cout << endl;

  istringstream sin("1.414 2.718 3.141 abc def, ghi");
  Point<3> p3, p4;
  sin >> p3; sin >> p4;
  cout << p3 << endl; cout << p4 << endl; 
  cout << endl;

  Point<3> p5{1, 1, 1};
  p5 *= 3; cout << p5 << endl;
  p5 /= 2; cout << p5 << endl;
  try {
    p5 /= 0; 
  } catch (exception &e) {
    cout << e.what() << endl;
  }
  cout << endl;

  Point<3> p6({2, 2, 2}), ones{1,1,1};
  p6 += ones; cout << p6 << ", " << ones << endl;
  p6 -= ones; cout << p6 << ", " << ones << endl;
  cout << endl;

  Point<3> twos{2,2,2};
  Point<3> p7 = twos + ones; Point<3> p8 = twos - ones;
  cout << p7 << ", " << p8 << endl;
  cout << endl;

  Point<3> p9 = 2 * twos; Point<3> p10 = twos * 3; 
  cout << p9 << ", " << p10 << endl;

  Point<3> p11 = twos / 3; Point<3> p12 = 2 * (twos - ones) / 3;
  cout << p11 << ", " << p12 << endl;


  return 0;
}
