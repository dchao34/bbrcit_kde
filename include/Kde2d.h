#ifndef BBRCITKDE_KDE2D_H__
#define BBRCITKDE_KDE2D_H__

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

  protected:

    double h1 = 1.0; double h2 = 1.0;
    std::vector<std::pair<double, double>> sample;

  private:

    double gauss2d(double, double, double, double);

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

#endif
