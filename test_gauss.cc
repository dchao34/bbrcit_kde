#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>

#include "gauss_legendre.h"
#include "Kde2d.h"
#include "ProdKde2d.h"

using namespace std;

double gauss2d(double x1, double x2, void *data) {
  return 1 / (2 * M_PI ) * exp(-0.5 * (x1 * x1 + x2 * x2));
}

double kde_wrapper(double x1, double x2, void *kde_obj_addr) {
  Kde2d *kde = (Kde2d *) kde_obj_addr;
  return (*kde)(x1, x2);
}

int main() {
  cout << gauss_legendre_2D_cube(10, gauss2d, nullptr, -2, 2, -2, 2) << endl;

  /*
  Kde2d kde("evttype1.dat");
  //kde.set_h1(0.0001);
  //kde.set_h2(0.001);
  kde.set_h1(0.001);
  kde.set_h2(0.01);

  //cout << gauss_legendre_2D_cube(128, kde_wrapper, &kde, -0.01, 0.1, 0.0, 1.0) << endl;
  cout << gauss_legendre_2D_cube(128, kde_wrapper, &kde, -0.02, 0.12, -0.1, 1.1) << endl;
  //cout << gauss_legendre_2D_cube(1024, kde_wrapper, &kde, -0.01, 0.1, 0.0, 1.0) << endl;
  //
  */

  ProdKde2d kde("evttype1.dat");
  kde.set_h1(0.001);
  kde.set_h2(0.01);
  cout << gauss_legendre_2D_cube(128, kde_wrapper, &kde, -0.02, 0.12, -0.1, 1.1) << endl;
  
}
