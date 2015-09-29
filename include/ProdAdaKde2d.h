#ifndef BBRCITKDE_PRODADAKDE2D_H__
#define BBRCITKDE_PRODADAKDE2D_H__

#define _USE_MATH_DEFINES
#include <cmath>

#include <string>
#include <vector>
#include <utility>

#include <ProdKde2d.h>

class ProdAdaKde2d : public ProdKde2d {

  public:

    // name of the file containing the data sample. 
    ProdAdaKde2d(std::string fname, 
                 double bw1, double bw2, 
                 double a1=0.5, double a2=0.5, 
                 unsigned r=20);

    ProdAdaKde2d() = default;
    ~ProdAdaKde2d() override = default;

    // evaluate the density estimate at a new point
    double operator()(double x1, double x2) override;
    double f1(double x) override { return evaluate_marginal(x, true); }
    double f2(double x) override { return evaluate_marginal(x, false); }

    // cross validation
    void cv(std::vector<double> &results, double h, bool cv_x1=true);

  private:

    double alpha1, alpha2;
    std::vector<std::pair<double,double>> bwidths;

    void compute_adaptive_bandwidths(
        std::vector<double> &results, bool dim1, unsigned r);

    double evaluate_marginal(double x, bool dim1);

};

inline double get_first_a(const std::pair<double,double> &p) { return p.first; }
inline double get_second_a(const std::pair<double,double> &p) { return p.second; }

#endif
