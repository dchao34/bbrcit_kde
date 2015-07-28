#ifndef __PRODKDE2D_H__
#define __PRODKDE2D_H__

#define _USE_MATH_DEFINES
#include <cmath>

#include <string>

#include "Kde2d.h"

class ProdKde2d : public Kde2d {

  public:

    // name of the file containing the data sample. 
    ProdKde2d(std::string fname) : Kde2d(fname) { compute_sample_stats(); };

    ProdKde2d() = default;
    ~ProdKde2d() override = default;

    // evaluate the density estimate at a new point
    double operator()(double x1, double x2) override;
    double f1(double x) { return evaluate_marginal(x, true); }
    double f2(double x) { return evaluate_marginal(x, false); }

    // compute silverman's bandwidth
    double silverman_h1() const { return pow(sample.size(), -1/6.0) * shat1; }
    double silverman_h2() const { return pow(sample.size(), -1/6.0) * shat2; }

    // perform cross validation
    void cv(std::ostream &os, const std::vector<double> &candidates, bool cv_h1=true);
    void cv(std::vector<double> &results, double h, bool cv_x1=true);

  private:

    double mhat1 = 0; double mhat2 = 0;
    double shat1 = 0; double shat2 = 0;
    
    double gauss_kernel_1d(double x, double s=1);
    double gauss_kernel_star(double x);
    void compute_sample_stats();

    double evaluate_marginal(double x, bool dim1);

};

// evaluates 1d gaussian kernel at x. (mean 0, stdev s)
inline double ProdKde2d::gauss_kernel_1d(double x, double s) {
  return 1/(s*sqrt(2*M_PI))*std::exp(-0.5*x*x/(s*s));
}

// evaluates the K* function. see Wasserman. 
inline double ProdKde2d::gauss_kernel_star(double x) {
  return gauss_kernel_1d(x, sqrt(2)) - 2 * gauss_kernel_1d(x, 1);
}

#endif
