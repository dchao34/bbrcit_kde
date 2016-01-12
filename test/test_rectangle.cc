#include <iostream>
#include <utility>
#include <Rectangle.h>

using namespace std;
using bbrcit::Rectangle;

int main() {

  Rectangle<2,double> r1, r2({0.0, 0.0}, {2.0, 2.0});
  cout << r1 << ", " << r2 << endl;
  swap(r1, r2);
  cout << r1 << ", " << r2 << endl;
  cout << endl;

  Rectangle<2,double> r3({1.0,1.0}, {3.0,3.0}); 
  Rectangle<2,double> r4({3.0, 3.0}, {4.0, 4.0});
  Rectangle<2,double> r5({0.0, 0.0}, {2.0, 2.0});
  cout << r3 << " " << r4 << " " << r5 << endl;
  cout << intersect(r5, r3) << " " << intersect(r5, r4) << endl;


  return 0;
}
