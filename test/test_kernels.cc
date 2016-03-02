#include <iostream>
#include <fstream>
#include <vector>

#include <Kernels/EpanechnikovKernel.h>
#include <Kernels/GaussianKernel.h>
#include <Kernels/EpanechnikovProductKernel.h>
#include <Kernels/GaussianProductKernel.h>
#include <DecoratedPoint.h>
#include <Attributes/KdeAttributes.h>

#include "kde_test_utils.h"

using namespace std;
using FloatType = double;
using bbrcit::EpanechnikovKernel;
using bbrcit::GaussianKernel;
using bbrcit::EpanechnikovProductKernel;
using bbrcit::GaussianProductKernel;
using bbrcit::KdeAttributes;
using Point1d = bbrcit::DecoratedPoint<1, KdeAttributes<FloatType>>;
using Point2d = bbrcit::DecoratedPoint<2, KdeAttributes<FloatType>>;

int main() {

  ofstream fout;
  vector<Point1d> grid1d; generate_1dgrid(grid1d, -5, 5, 1000);
  vector<Point2d> grid2d; generate_2dgrid(grid2d, -5, 5, 1000, -5, 5, 1000);

  // test: EpanechnikovKernel
  EpanechnikovKernel<1,FloatType> ep1d, ep1d_narrow(0.5), ep1d_wide(2);
  cout << ep1d.dim() << " " << ep1d.normalization() << " (c.f. 1, 0.75) " << std::endl;
  cout << ep1d_narrow.dim() << " " << ep1d_narrow.normalization() << " (c.f. 1, 1.5) " << std::endl;
  cout << ep1d_wide.dim() << " " << ep1d_wide.normalization() << " (c.f. 1, 0.375) " << std::endl;
  fout.open("test_epanechnikov1d.csv");
  for (const auto &p : grid1d) { 
    fout << p[0] << " " << ep1d.eval(p);
    fout << " " << ep1d_narrow.eval(p);
    fout << " " << ep1d_wide.eval(p);
    fout << endl; 
  }
  fout.close();

  EpanechnikovKernel<2,FloatType> ep2d;
  cout << ep2d.dim() << " " << ep2d.normalization();
  cout << " (c.f. 2, "<< 0.5 * (1/M_PI) * 4 << ") " << std::endl;
  fout.open("test_epanechnikov2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = ep2d.eval(p); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  // test: GaussianKernel
  GaussianKernel<1,FloatType> gauss1d, gauss1d_narrow(0.5), gauss1d_wide(2);
  cout << gauss1d.dim() << " " << gauss1d.normalization();
  cout << " (c.f. 1, "<< std::pow(2*M_PI, -0.5)<< ") " << std::endl;
  cout << gauss1d_narrow.dim() << " " << gauss1d_narrow.normalization();
  cout << " (c.f. 1, "<< std::pow(2*M_PI, -0.5) / 0.5 << ") " << std::endl;
  cout << gauss1d_wide.dim() << " " << gauss1d_wide.normalization();
  cout << " (c.f. 1, "<< std::pow(2*M_PI, -0.5) / 2.0 << ") " << std::endl;
  fout.open("test_gauss1d.csv");
  for (const auto &p : grid1d) { 
    fout << p[0] << " " << gauss1d.eval(p);
    fout << " " << gauss1d_narrow.eval(p);
    fout << " " << gauss1d_wide.eval(p);
    fout << endl; 
  }
  fout.close();

  GaussianKernel<2,FloatType> gauss2d;
  cout << gauss2d.dim() << " " << gauss2d.normalization();
  cout << " (c.f. 2, "<< std::pow(2*M_PI, -1.0)<< ") " << std::endl;
  fout.open("test_gauss2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = gauss2d.eval(p); 
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

  GaussianProductKernel<2> gaussprod2d({1.2, 1.0});
  cout << gaussprod2d.dim() << " " << gaussprod2d.normalization();
  cout << " (c.f. 2, "<< std::pow(2*M_PI, -1.0) / (1.2*1.0) << ") " << std::endl;
  gaussprod2d.set_bandwidth(0, 1.5);
  cout << gaussprod2d.dim() << " " << gaussprod2d.normalization();
  cout << " (c.f. 2, "<< std::pow(2*M_PI, -1.0) / (1.5*1.0) << ") " << std::endl;
  gaussprod2d.set_bandwidths({1.2, 0.7});
  cout << gaussprod2d.dim() << " " << gaussprod2d.normalization();
  cout << " (c.f. 2, "<< std::pow(2*M_PI, -1.0) / (1.2*0.7) << ") " << std::endl;
  fout.open("test_gaussprod2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = gaussprod2d.eval(p); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  // test: EpanechnikovProductKernel
  EpanechnikovKernel<1> epprod1d;
  cout << epprod1d.dim() << " " << epprod1d.normalization();
  cout << " (c.f. 1, "<< std::pow(0.75, 1)<< ") " << std::endl;
  fout.open("test_epprod1d.csv");
  for (const auto &p : grid1d) { fout << p[0] << " " << epprod1d.eval(p) << endl; }
  fout.close();

  EpanechnikovProductKernel<2> epprod2d({1.2, 1.0});
  cout << epprod2d.dim() << " " << epprod2d.normalization();
  cout << " (c.f. 2, "<< std::pow(0.75, 2) / (1.2*1.0) << ") " << std::endl;
  epprod2d.set_bandwidth(0, 1.5);
  cout << epprod2d.dim() << " " << epprod2d.normalization();
  cout << " (c.f. 2, "<< std::pow(0.75, 2) / (1.5*1.0) << ") " << std::endl;
  epprod2d.set_bandwidths({1.2, 0.7});
  cout << epprod2d.dim() << " " << epprod2d.normalization();
  cout << " (c.f. 2, "<< std::pow(0.75, 2) / (1.2*0.7) << ") " << std::endl;
  fout.open("test_epprod2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = epprod2d.eval(p); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  return 0;
}
