#ifndef BBRCITKDE_PRODKDE2D_H__
#define BBRCITKDE_PRODKDE2D_H__

#define _USE_MATH_DEFINES
#include <cmath>

#include <string>
#include <vector>
#include <utility>
#include <complex>

#include <Kde2d.h>

class ProdKde2d : public Kde2d {

  protected:
    using vp_it = std::vector<std::pair<double,double>>::iterator;
    using cplex = std::complex<double>;
    using vc_size_t = std::vector<cplex>::size_type;
    using vp_size_t = std::vector<std::pair<double,double>>::size_type;

  public:

    // name of the file containing the data sample. 
    ProdKde2d(std::string fname, double bw1=1.0, double bw2=1.0) 
      : Kde2d(fname, bw1, bw2) { compute_sample_stats(); };

    ProdKde2d() = default;
    virtual ~ProdKde2d() override = default;

    // evaluate the density estimate at a new point
    virtual double operator()(double x1, double x2) override;
    virtual double f1(double x) { return evaluate_marginal(x, true); }
    virtual double f2(double x) { return evaluate_marginal(x, false); }
    double f1_se(double x) { return evaluate_marginal_se(x, true); }
    double f2_se(double x) { return evaluate_marginal_se(x, false); }
    void grid_evaluate_marginal(std::vector<std::pair<double,double>> &results, 
                                bool dim1, unsigned r=20);

    // compute silverman's bandwidth
    double silverman_h1() const { return pow(sample.size(), -1/6.0) * shat1; }
    double silverman_h2() const { return pow(sample.size(), -1/6.0) * shat2; }

    // perform cross validation
    void cv(std::vector<double> &results, double h, bool cv_x1=true);
    void fcv(std::vector<double> &results, double h, unsigned r, bool cv_x1=true);


  protected:
    double gauss_kernel_1d(double x, double s=1);
    double gauss_kernel_star(double x);

  private:

    double evaluate_marginal(double x, bool dim1);
    double evaluate_marginal_se(double x, bool dim1);

    // silverman bandwidth helpers
    double mhat1 = 0; double mhat2 = 0;
    double shat1 = 0; double shat2 = 0;
    void compute_sample_stats();

    // fft helpers
    void deduce_fft_grid_constants(double&, double&, double&, vc_size_t, double, bool);
    void discretize_data(std::vector<double>&, double, double, double,
                         vc_size_t, vc_size_t, bool);
    void fft_forward(const std::vector<double>&, std::vector<cplex>&, unsigned);

};

// evaluates 1d gaussian kernel at x. (mean 0, stdev s)
inline double ProdKde2d::gauss_kernel_1d(double x, double s) {
  return 1/(s*sqrt(2*M_PI))*std::exp(-0.5*x*x/(s*s));
}

// evaluates the K* function. see Wasserman. 
inline double ProdKde2d::gauss_kernel_star(double x) {
  return gauss_kernel_1d(x, sqrt(2)) - 2 * gauss_kernel_1d(x, 1);
}

inline double get_first(std::vector<std::pair<double,double>>::iterator p) { return p->first; }
inline double get_second(std::vector<std::pair<double,double>>::iterator p) { return p->first; }

#endif
