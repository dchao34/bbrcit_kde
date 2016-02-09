#include <iostream>
#include <fstream>
#include <vector>

#include <EpanechnikovKernel.h>
#include <GaussianKernel.h>
#include <GaussianProductKernel.h>
#include <DecoratedPoint.h>
#include <Attributes/KdeAttributes.h>

#include "kde_test_utils.h"

using namespace std;
using bbrcit::EpanechnikovKernel;
using bbrcit::GaussianKernel;
using bbrcit::GaussianProductKernel;
using bbrcit::KdeAttributes;
using Point1d = bbrcit::DecoratedPoint<1, KdeAttributes<double>>;
using Point2d = bbrcit::DecoratedPoint<2, KdeAttributes<double>>;

int main() {

  ofstream fout;
  vector<Point1d> grid1d; generate_1dgrid(grid1d, -5, 5, 1000);
  vector<Point2d> grid2d; generate_2dgrid(grid2d, -5, 5, 1000, -5, 5, 1000);

  // test: EpanechnikovKernel
  EpanechnikovKernel<1> ep1d;
  cout << ep1d.dim() << " " << ep1d.normalization();
  cout << " (c.f. 1, 0.75) " << std::endl;
  fout.open("test_epanechnikov1d.csv");
  for (const auto &p : grid1d) { fout << p[0] << " " << ep1d.eval(p) << endl; }
  fout.close();

  EpanechnikovKernel<2> ep2d;
  cout << ep2d.dim() << " " << ep2d.normalization();
  cout << " (c.f. 2, "<< 0.5 * (1/M_PI) * 4 << ") " << std::endl;
  fout.open("test_epanechnikov2d.csv");
  for (auto &p : grid2d) { 
    double result = ep2d.eval(p); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  // test: GaussianKernel
  GaussianKernel<1> gauss1d;
  cout << gauss1d.dim() << " " << gauss1d.normalization();
  cout << " (c.f. 1, "<< std::pow(2*M_PI, -0.5)<< ") " << std::endl;
  fout.open("test_gauss1d.csv");
  for (const auto &p : grid1d) { fout << p[0] << " " << gauss1d.eval(p) << endl; }
  fout.close();

  GaussianKernel<2> gauss2d;
  cout << gauss2d.dim() << " " << gauss2d.normalization();
  cout << " (c.f. 2, "<< std::pow(2*M_PI, -1.0)<< ") " << std::endl;
  fout.open("test_gauss2d.csv");
  for (auto &p : grid2d) { 
    double result = gauss2d.eval(p); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  // test: GaussianProductKernel
  GaussianProductKernel<1> gaussprod1d;
  cout << gaussprod1d.dim() << " " << gaussprod1d.normalization();
  cout << " (c.f. 1, "<< std::pow(2*M_PI, -0.5)<< ") " << std::endl;
  fout.open("test_gaussprod1d.csv");
  for (const auto &p : grid1d) { fout << p[0] << " " << gaussprod1d.eval(p) << endl; }
  fout.close();

  GaussianProductKernel<2> gaussprod2d({1.2, 0.7});
  cout << gaussprod2d.dim() << " " << gaussprod2d.normalization();
  cout << " (c.f. 2, "<< std::pow(2*M_PI, -1.0) / (1.2*0.7) << ") " << std::endl;
  fout.open("test_gaussprod2d.csv");
  for (auto &p : grid2d) { 
    double result = gaussprod2d.eval(p); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  return 0;
}
