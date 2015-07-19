#include <iostream>
#include <fstream>
#include <sstream>

#include "Kde2d.h"
#include "ProdKde2d.h"

using namespace std;

int main() {

  ProdKde2d kde("evttype3.dat");

  double h1 = 0.01, h2 = 0.1;
  kde.set_h1(h1); kde.set_h2(h2);

  cout << "start" << endl;
  cout << kde(0.1, 0.01) << endl;
  cout << "end" << endl;

  return 0;
}
