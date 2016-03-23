#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include <Kernels/EpanechnikovKernel.h>
#include <Kernels/GaussianKernel.h>
#include <Kernels/EpanechnikovProductKernel2d.h>
#include <Kernels/GaussianProductKernel2d.h>
#include <Kernels/EpanechnikovConvKernel1d.h>
#include <Kernels/EpanechnikovProductConvKernel2d.h>
#include <Kernels/GaussianConvKernel1d.h>
#include <Kernels/GaussianProductConvKernel2d.h>

#include <DecoratedPoint.h>
#include <Attributes/AdaKdeAttributes.h>

#include "kde_test_utils.h"

using namespace std;
using FloatType = float;

using bbrcit::EpanechnikovKernel;
using bbrcit::GaussianKernel;
using bbrcit::EpanechnikovProductKernel2d;
using bbrcit::GaussianProductKernel2d;
using bbrcit::EpanechnikovConvKernel1d;
using bbrcit::GaussianConvKernel1d;
using bbrcit::EpanechnikovProductConvKernel2d;
using bbrcit::GaussianProductConvKernel2d;

using bbrcit::AdaKdeAttributes;
using PointType1d = bbrcit::DecoratedPoint<1, AdaKdeAttributes<FloatType>>;
using PointType2d = bbrcit::DecoratedPoint<2, AdaKdeAttributes<FloatType>>;

const PointType1d origin1d; 
const PointType2d origin2d;

int main() {

  std::mt19937 e;
  vector<FloatType> sim_p;

  ofstream fout;
  vector<PointType1d> grid1d; generate_1dgrid(grid1d, -5, 5, 1000);
  vector<PointType2d> grid2d; generate_2dgrid(grid2d, -5, 5, 1000, -5, 5, 1000);

  // test: EpanechnikovKernel
  EpanechnikovKernel<1,FloatType> ep1d, ep1d_narrow(0.5), ep1d_wide(2);
  cout << ep1d.dim() << " " << ep1d.normalization() << " (c.f. 1, 0.75) " << std::endl;
  cout << ep1d_narrow.dim() << " " << ep1d_narrow.normalization() << " (c.f. 1, 1.5) " << std::endl;
  cout << ep1d_wide.dim() << " " << ep1d_wide.normalization() << " (c.f. 1, 0.375) " << std::endl;

  fout.open("test_epanechnikov1d.csv");
  for (const auto &p : grid1d) { 
    fout << p[0] << " " << ep1d.normalization() * ep1d.unnormalized_eval(p, origin1d);
    fout << " " << ep1d_narrow.normalization() * ep1d_narrow.unnormalized_eval(p, origin1d);
    fout << " " << ep1d_wide.normalization() * ep1d_wide.unnormalized_eval(p, origin1d);
    fout << endl; 
  }
  fout.close();

  EpanechnikovKernel<2,FloatType> ep2d;
  cout << ep2d.dim() << " " << ep2d.normalization();
  cout << " (c.f. 2, "<< 0.5 * (1/M_PI) * 4 << ") " << std::endl;
  fout.open("test_epanechnikov2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = ep2d.normalization() * ep2d.unnormalized_eval(p, origin2d); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  // simulate
  fout.open("test_epanechnikov1d_sim.csv");
  for (size_t i = 0; i < 100000; ++i) {
    ep1d.simulate(e, sim_p);
    fout << sim_p[0] << " ";
    ep1d.simulate(e, sim_p, 0.5);
    fout << sim_p[0] << " ";
    ep1d.simulate(e, sim_p, 2.0);
    fout << sim_p[0] << endl;
  }
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
    fout << p[0] << " " << gauss1d.normalization() * gauss1d.unnormalized_eval(p, origin1d);
    fout << " " << gauss1d_narrow.normalization() * gauss1d_narrow.unnormalized_eval(p, origin1d);
    fout << " " << gauss1d_wide.normalization() * gauss1d_wide.unnormalized_eval(p, origin1d);
    fout << endl; 
  }
  fout.close();

  // simulate
  fout.open("test_gauss1d_sim.csv");
  for (size_t i = 0; i < 100000; ++i) {
    gauss1d.simulate(e, sim_p);
    fout << sim_p[0] << " ";
    gauss1d.simulate(e, sim_p, 0.5);
    fout << sim_p[0] << " ";
    gauss1d.simulate(e, sim_p, 2.0);
    fout << sim_p[0] << endl;
  }
  fout.close();

  GaussianKernel<2,FloatType> gauss2d;
  cout << gauss2d.dim() << " " << gauss2d.normalization();
  cout << " (c.f. 2, "<< std::pow(2*M_PI, -1.0)<< ") " << std::endl;
  fout.open("test_gauss2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = gauss2d.normalization() * gauss2d.unnormalized_eval(p, origin2d); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  // simulate
  fout.open("test_gauss2d_sim.csv");
  for (size_t i = 0; i < 100000; ++i) {
    gauss2d.simulate(e, sim_p);
    fout << sim_p[0] << " " << sim_p[1] << endl;
  }
  fout.close();

  // test: GaussianProductKernel
  GaussianProductKernel2d<FloatType> gaussprod2d(1.2, 1.0);
  cout << gaussprod2d.normalization();
  cout << " (c.f. "<< std::pow(2*M_PI, -1.0) / (1.2*1.0) << ") " << std::endl;
  gaussprod2d.set_hx(1.5);
  cout << gaussprod2d.normalization();
  cout << " (c.f. "<< std::pow(2*M_PI, -1.0) / (1.5*1.0) << ") " << std::endl;
  gaussprod2d.set_bandwidths(1.2, 0.7);
  cout << gaussprod2d.normalization();
  cout << " (c.f. "<< std::pow(2*M_PI, -1.0) / (1.2*0.7) << ") " << std::endl;
  fout.open("test_gaussprod2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = gaussprod2d.normalization() * gaussprod2d.unnormalized_eval(p, origin2d); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  // simulate
  fout.open("test_gaussprod2d_sim.csv");
  for (size_t i = 0; i < 100000; ++i) {
    gaussprod2d.simulate(e, sim_p);
    fout << sim_p[0] << " " << sim_p[1] << endl;
  }
  fout.close();

  // test: EpanechnikovProductKernel2d
  EpanechnikovProductKernel2d<FloatType> epprod2d(1.2, 1.0);
  cout << epprod2d.normalization();
  cout << " (c.f. "<< std::pow(0.75, 2) / (1.2*1.0) << ") " << std::endl;
  epprod2d.set_hx(1.5);
  cout << epprod2d.normalization();
  cout << " (c.f. "<< std::pow(0.75, 2) / (1.5*1.0) << ") " << std::endl;
  epprod2d.set_bandwidths(1.2, 0.7);
  cout << epprod2d.normalization();
  cout << " (c.f. "<< std::pow(0.75, 2) / (1.2*0.7) << ") " << std::endl;
  fout.open("test_epprod2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = epprod2d.normalization() * epprod2d.unnormalized_eval(p, origin2d); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  // simulate
  fout.open("test_epprod2d_sim.csv");
  for (size_t i = 0; i < 100000; ++i) {
    epprod2d.simulate(e, sim_p);
    fout << sim_p[0] << " " << sim_p[1] << endl;
  }
  fout.close();

  // test: EpanechnikovConvKernel1d
  EpanechnikovConvKernel1d<FloatType> epconv1d, epconv1d_narrow(0.5), epconv1d_wide(2.0);

  fout.open("test_epconv1d.csv");
  for (const auto &p : grid1d) { 
    fout << p[0] << " " << epconv1d.normalization() * epconv1d.unnormalized_eval(p, origin1d);
    fout << " " << epconv1d_narrow.normalization() * epconv1d_narrow.unnormalized_eval(p, origin1d);
    fout << " " << epconv1d_wide.normalization() * epconv1d_wide.unnormalized_eval(p, origin1d);
    fout << endl; 
  }
  fout.close();

  // test: EpanechnikovProductConvKernel2d
  EpanechnikovProductConvKernel2d<FloatType> epconv2d(2.0, 1.0);
  fout.open("test_epconv2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = epconv2d.normalization() * epconv2d.unnormalized_eval(p, origin2d); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  // test: GaussianConvKernel1d
  GaussianConvKernel1d<FloatType> gaussconv1d, gaussconv1d_narrow(0.5), gaussconv1d_wide(2.0);

  fout.open("test_gaussconv1d.csv");
  for (const auto &p : grid1d) { 
    fout << p[0] << " " << gaussconv1d.normalization() * gaussconv1d.unnormalized_eval(p, origin1d);
    fout << " " << gaussconv1d_narrow.normalization() * gaussconv1d_narrow.unnormalized_eval(p, origin1d);
    fout << " " << gaussconv1d_wide.normalization() * gaussconv1d_wide.unnormalized_eval(p, origin1d);
    fout << endl; 
  }
  fout.close();

  // test: GaussianProductConvKernel2d
  GaussianProductConvKernel2d<FloatType> gaussconv2d(1.2, 1.0);
  fout.open("test_gaussconv2d.csv");
  for (auto &p : grid2d) { 
    FloatType result = gaussconv2d.normalization() * gaussconv2d.unnormalized_eval(p, origin2d); 
    p.attributes().set_lower(result);
    p.attributes().set_upper(result);
  };
  write_kde2d_result(fout, grid2d, -5, 5, 1000, -5, 5, 1000);
  fout.close();

  return 0;
}
