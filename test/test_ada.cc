#include "ProdAdaKde2d.h"

using namespace std; 

int main() {
  ProdAdaKde2d kde("data/evttype1.dither.csv", 1.97e-3, 1.95e-2, 0.5, 0.5);
  cout << kde(0.03, 0.3) << endl;
  cout << kde.f1(0.03) << endl;
  cout << kde.f2(0.3) << endl;
  return 0;
}
