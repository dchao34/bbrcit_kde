#ifndef __KDE2D_H__
#define __KDE2D_H__

#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <string>
#include <vector>
#include <utility>

class Kde2d {

  protected:

    using sample_no = std::vector<std::pair<double, double>>::size_type;

  public:

    // name of the file containing the data sample. 
    Kde2d(std::string, double bw1 = 1.0, double bw2 = 1.0);
    Kde2d() = default;
    virtual ~Kde2d() = default;

    // set smoothing parameters
    void set_h1(double h) { h1 = h; }
    void set_h2(double h) { h2 = h; }

    // evaluate the density estimate at a new point
    virtual double operator()(double x1, double x2);

    // perform cross validation
    void cv(std::ostream &os, 
            const std::vector<std::pair<double, double>> &candidates, 
            double x_low, double x_high, double y_low, double y_high,
            int qgauss_n);

  protected:

    double h1 = 1.0; double h2 = 1.0;
    std::vector<std::pair<double, double>> sample;

  private:

    double gauss2d(double, double, double, double);
    double excluded_eval(double x1, double x2, sample_no bidx, sample_no eidx);
    double eval2(double x1, double x2);

    friend double f2(double x1, double x2, void *kde_obj_addr);

};

// evaluates the 2d gaussian density at (x1, x2) with mean = (m1, m2)
// and covariance ((h1*h1, 0), (0, h2*h2));
inline double Kde2d::gauss2d(
    double x1, double x2, 
    double m1, double m2) { 
  return 1 / (2 * M_PI * h1 * h2) *
        exp(
            -0.5 * (
              (x1-m1) * (x1-m1) / (h1 * h1) + 
              (x2-m2) * (x2-m2) / (h2 * h2)
            )
        );
}

// wrapper to conform to the API required by the gaussian quadrature 
// code written by Pavel Holoborodko:
// http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/
double f2(double x1, double x2, void *kde_obj_addr);

// evaluate the square of the  estimate at a new point
inline double Kde2d::eval2(double x1, double x2) { 
  double res = (*this)(x1, x2); 
  return res * res; 
}

#endif
