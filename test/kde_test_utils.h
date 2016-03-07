#ifndef BBRCITKDE_KDETESTUTILS_H__
#define BBRCITKDE_KDETESTUTILS_H__

#define _USE_MATH_DEFINES

#ifdef __CUDACC__
#define CUDA_UNUSED __attribute__((unused))
#else 
#define CUDA_UNUSED
#endif

#include <vector>
#include <cmath>
#include <ostream>
#include <random>
#include <algorithm>

void transform(double &x, double &y,
               double mx, double my,
               double sx, double sy,
               double degrees) {
  x *= sx; y *= sy;
  x += mx; y += my;
  double costh = std::cos(degrees * M_PI / 180.);
  double sinth = std::sin(degrees * M_PI / 180.);
  x = costh*x - sinth*y;
  y = sinth*x + costh*y;
}

template <typename RND_E, typename PointT>
void generate_bimodal_gaussian(
    RND_E &e, std::vector<PointT> &data, int n_samples, 
    double mx1, double sx1, double mx2, double sx2,
    double p1 = 0.5, double p1_wgt = 1.0) {

  data.clear(); data.reserve(n_samples);

  std::normal_distribution<> gaussian(0.0, 1.0);
  std::uniform_real_distribution<> uniform(0.0, 1.0);

  for (int i = 0; i < n_samples; ++i) {
    double x = gaussian(e), u = uniform(e);
    double CUDA_UNUSED w = 1.0;
    if (u < p1) {
      x *= sx1; x += mx1;
      w = p1_wgt;
    } else {
      x *= sx2; x += mx2;
    }
    data.push_back({{x}, {w}});
  }
}

template <typename RND_E, typename PointT>
void generate_bimodal_gaussian(
    RND_E &e, std::vector<PointT> &data, int n_samples, 
    double mx1, double my1, double sx1, double sy1, double dgr1,
    double mx2, double my2, double sx2, double sy2, double dgr2, 
    double p1 = 0.5, double p1_wgt = 1.0) {

  data.clear(); data.reserve(n_samples);

  std::normal_distribution<> gaussian(0.0, 1.0);
  std::uniform_real_distribution<> uniform(0.0, 1.0);

  for (int i = 0; i < n_samples; ++i) {
    double x = gaussian(e), y = gaussian(e), u = uniform(e);
    double CUDA_UNUSED w = 1.0;
    if (u < p1) {
      transform(x, y, mx1, my1, sx1, sy1, dgr1);
      w = p1_wgt;
    } else {
      transform(x, y, mx2, my2, sx2, sy2, dgr2);
    }
    data.push_back({{x,y}, {w}});
  }
}

template <typename PointT>
void generate_1dgrid(
    std::vector<PointT> &grid, 
    double start_x, double end_x, double steps_x) {
  
  grid.clear();

  double delta_x = (end_x - start_x) / steps_x;

  double x_coord;
  for (int i = 0; i < steps_x; ++i) {
    x_coord = start_x + i * delta_x;
    PointT p = {{ x_coord }};
    grid.push_back(p);
  }

}

template <typename PointT>
void generate_2dgrid(
    std::vector<PointT> &grid, 
    double start_x, double end_x, double steps_x,
    double start_y, double end_y, double steps_y) {
  
  grid.clear();

  double delta_x = (end_x - start_x) / steps_x;
  double delta_y = (end_y - start_y) / steps_y;

  double CUDA_UNUSED x_coord, y_coord;
  for (int j = 0; j < steps_y; ++j) {

    y_coord = start_y + j * delta_y;

    for (int i = 0; i < steps_x; ++i) {

      x_coord = start_x + i * delta_x;

      grid.push_back({{x_coord, y_coord}});

    }
  }


}

template<typename PointT>
bool ReverseExactLexicoLess(const PointT &lhs, const PointT &rhs) {
  int i = 0; while (i < lhs.dim() && lhs[lhs.dim()-i-1] == rhs[lhs.dim()-i-1]) { ++i; }
  return i != lhs.dim() && lhs[lhs.dim()-i-1] < rhs[lhs.dim()-i-1];
}

template<typename PointT>
void write_scatter_data(std::ostream &os, const std::vector<PointT> &data) {
  for (const auto &p : data) { os << p[0] << " " << p[1] << std::endl; }
}

template<typename PointT>
void write_kde2d_result(
    std::ostream &os, std::vector<PointT> &queries,
    double start_x, double end_x, int steps_x, 
    double start_y, double end_y, int steps_y) {

  double delta_x = (end_x - start_x) / steps_x;
  double delta_y = (end_y - start_y) / steps_y;
  for (int i = 0; i < steps_x; ++i) { os << start_x + i * delta_x << " "; } 
  os << std::endl;
  for (int i = 0; i < steps_y; ++i) { os << start_y + i * delta_y << " "; } 
  os << std::endl;

  std::sort(queries.begin(), queries.end(), 
            ReverseExactLexicoLess<PointT>);
  for (size_t i = 0; i < queries.size(); ++i) {
    if (i % steps_x == 0 && i) { os << std::endl; }
    os << queries[i].attributes().middle() << " ";
  }
  os << std::endl;
}




#endif
