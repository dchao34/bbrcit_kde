#include <iostream>
#include <vector>
#include <Point.h>

using namespace std;

int main() {

  char status_c;
  cout << "Point1d Tests" << endl;
  cout << "-------------" << endl;
  cout << endl;

  cout << "Test 1: Constructors. " << endl;

  // 1.1 default constructor
  Point1d p1;  // p1 = 0.0
  status_c = (p1.x() == 0) ? '.' : 'x';
  cout << status_c << flush;

  // 1.2 double constructor
  Point1d p2(3.0); // p2 = 3.0
  status_c = (p2.x() == 3.0) ? '.' : 'x';
  cout << status_c << flush;

  cout << endl;

  cout << "Test 2: Copy-control. " << endl;

  // 2.1 copy constructor
  Point1d p3(p2); // p2 = 3.0, p3 = 3.0
  status_c = (p3.x() == 3.0) ? '.' : 'x';
  cout << status_c << flush;

  // 2.2 copy assignment
  p3 = p1; // p1 = 0.0, p3 = 0.0
  status_c = (p3.x() == 0.0) ? '.' : 'x';
  cout << status_c << flush;

  // 2.3 implicit conversion
  p1 = 4.0; // p1 = 4.0
  status_c = (p1.x() == 4.0) ? '.' : 'x';
  cout << status_c << flush;
  
  // 2.4 self assignment
  p2 = p2; // p2 = 3.0
  status_c = (p2.x() == 3.0) ? '.' : 'x';
  cout << status_c << flush;

  // 2.5 chain assignment
  p2 = p3 = p1;  // p1 = 4.0, p2 = 4.0, p3 = 4.0
  if (p1.x() == 4.0 && p2.x() == 4.0 && p3.x() == 4.0) {
    status_c = '.';
  } else {
    status_c = 'x';
  }
  cout << status_c << flush;

  cout << endl;

  cout << "Test 3: Element access. " << endl;

  // 3.1 access
  const Point1d p4(p3); // p3 = 4.0, p4 = 4.0
  if (p4[0] == 4.0 && p3[0] == 4.0) { 
    status_c = '.';
  } else {
    status_c = 'x';
  }
  cout << status_c << flush;

  // 3.2 mutate
  ++p3[0]; // p3: 4.0 -> 5.0
  status_c = (p3.x() == 5.0) ? '.' : 'x';
  cout << status_c << flush;

  cout << endl;

  cout << "Test 4: Scaler Arithmetic. " << endl;

  // scalar multiplication
  // ---------------------

  // 4.1 operator*=
  p3[0] = 4.0;
  p3 *= 2; // p3 = 8.0
  status_c = (p3.x() == 8.0) ? '.' : 'x';
  cout << status_c << flush;

  // 4.2 operator*=
  p3 *= -2; // p3 = -16.0
  status_c = (p3.x() == -16.0) ? '.' : 'x';
  cout << status_c << flush;

  // 4.3 operator*
  Point1d p5 = p3 * -1.0;
  status_c = (p5.x() == 16.0) ? '.' : 'x';
  cout << status_c << flush;

  // 4.4 operator*
  p5 = 2.0 * p3;
  status_c = (p5.x() == -32.0) ? '.' : 'x';
  cout << status_c << flush;

  // 4.5 operator*
  p5 = -2.0 * p5;
  status_c = (p5.x() == 64.0) ? '.' : 'x';
  cout << status_c << flush;

  // scalar division
  // ---------------

  // 4.6 operator/=
  p5 /= 2.0;
  status_c = (p5.x() == 32.0) ? '.' : 'x';
  cout << status_c << flush;

  // 4.7 operator/=
  p5 /= -2.0;
  status_c = (p5.x() == -16.0) ? '.' : 'x';
  cout << status_c << flush;

  // 4.8 operator/
  p5 = p5 / -2.0;
  status_c = (p5.x() == 8.0) ? '.' : 'x';
  cout << status_c << flush;

  cout << endl;

  cout << "Test 5: Vector Arithmetic. " << endl;
  
  // 5.1 operator+=
  Point1d p6(1.0); p5[0] = 5.0;
  p5 += p6;
  if (p6.x() == 1.0 && p5.x() == 6.0) { 
    status_c = '.';
  } else {
    status_c = 'x';
  }
  cout << status_c << flush;
  
  // 5.2 operator-=
  p5 -= p6;
  if (p6.x() == 1.0 && p5.x() == 5.0) { 
    status_c = '.';
  } else {
    status_c = 'x';
  }
  cout << status_c << flush;
  
  // 5.3 operator+
  Point1d p7 = p5 + p6;
  if (p7.x() == 6.0 && p6.x() == 1.0 && p5.x() == 5.0) { 
    status_c = '.';
  } else {
    status_c = 'x';
  }
  cout << status_c << flush;
  
  // 5.3 operator-
  p7 = p5 - p6;
  if (p7.x() == 4.0 && p6.x() == 1.0 && p5.x() == 5.0) { 
    status_c = '.';
  } else {
    status_c = 'x';
  }
  cout << status_c << flush;

  cout << endl;


  cout << "Test 6: Comparison operators. " << endl;
  p1[0] = 1.0; p2[0] = 1.0; p3[0] = 2.0;


  // 5.1 operator==
  status_c = (p1 == p1 && p1 == p2 && !(p2 == p3) ) ? '.' : 'x';
  cout << status_c << flush;


  // 5.2 operator!=
  status_c = (!(p1 != p1) && !(p1 != p2) && p2 != p3) ? '.' : 'x';
  cout << status_c << flush;


  // 5.3 operator<
  status_c = (!(p1 < p1) && !(p3 < p2) && p2 < p3) ? '.' : 'x';
  cout << status_c << flush;


  // 5.4 operator<=
  status_c = (p1 <= p1 && !(p3 <= p2) && p2 <= p3) ? '.' : 'x';
  cout << status_c << flush;


  // 5.5 operator>
  status_c = (!(p1 > p1) && p3 > p2 && !(p2 > p3)) ? '.' : 'x';
  cout << status_c << flush;


  // 5.6 operator<=
  status_c = (p1 >= p1 && p3 >= p2 && !(p2 >= p3)) ? '.' : 'x';
  cout << status_c << flush;

  cout << endl;

  cout << "Test 7: Vector operations. " << endl;

  // 7.1 norm()
  p1[0] = -5.0;
  status_c = (p1.norm() == 5.0) ? '.' : 'x';
  cout << status_c << flush;

  // 7.2 norm2()
  status_c = (p1.norm2() == 25.0) ? '.' : 'x';
  cout << status_c << flush;

  cout << endl;

  cout << "Test 8: Container compliance. " << endl;

  // 7.1 sort in vector
  p1[0] = 1.0; p2[0] = 2.0; p3[0] = 3.0; // p4 already 4.0
  p5[0] = 5.0; p6[0] = 6.0;
  vector<Point1d> v;
  v.push_back(p2); v.push_back(p1);
  v.push_back(p4); v.push_back(p3);
  v.push_back(p6); v.push_back(p5);
  sort(v.begin(), v.end());
  int i; for (i = 1; i < v.size() && v[i][0] > v[i-1][0]; ++i) ;
  status_c = (i == v.size()) ? '.' : 'x';
  cout << status_c << flush;

  
  cout << endl;

  return 0;
}
