#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <Point.h>

using namespace std;

int main() {
  Point<3> p1;
  cout << p1 << endl;
  Point<3> p2 = { 1.0, 2.0, 3.2 };
  cout << p2 << endl;

  istringstream sin("1.414 2.718 3.141 abc def, ghi");
  Point<3> p3, p4;
  sin >> p3;
  cout << p3 << endl; 
  sin >> p4;
  cout << p4 << endl; 

  vector<double> v(3, 0);
  array<double, 3> a = {{ 0, 0, 0}};
  cout << sizeof(v) << " " << sizeof(a) << endl;

  return 0;
}
