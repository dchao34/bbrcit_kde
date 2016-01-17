#include <iostream>
#include <vector>
#include <random>

#include <Point.h>
#include <DecoratedPoint.h>
#include <FloatUtils.h>

using namespace std;

using bbrcit::Point;
using bbrcit::DecoratedPoint;
using Point2d = Point<2>;
using WPoint2d = DecoratedPoint<2>;

template<typename P2dT>
vector<P2dT> generate_grid(
    double start_x, double end_x, 
    double start_y, double end_y, 
    int N_x, int N_y) {
  vector<P2dT> result;
  double delta_x = (end_x - start_x) / N_x, delta_y = (end_y - start_y) / N_y;
  for (int i = 0; i < N_x; ++i) { 
    for (int j = 0; j < N_y; ++j) {
      result.push_back({{start_x + i*delta_x, start_y + j*delta_y}});
    }
  }
  return result;
}

int main() {

  cout << endl;

  vector<Point2d> points = generate_grid<Point2d>(0, 1, 0, 1, 3, 3);
  default_random_engine e;
  shuffle(points.begin(), points.end(), e);
  sort(points.begin(), points.end(), bbrcit::ExactLexicoLess<Point2d>);
  for (auto const &e : points) { cout << e << endl; }

  cout << endl;

  vector<WPoint2d> wpoints = generate_grid<WPoint2d>(0, 1, 0, 1, 3, 3);
  shuffle(wpoints.begin(), wpoints.end(), e);
  sort(wpoints.begin(), wpoints.end(), bbrcit::ExactLexicoLess<WPoint2d>);
  for (auto const &e : wpoints) { cout << e << endl; }

  cout << endl;

  WPoint2d wp0 = {{ 0.0, 1.0 } , 2.0 };
  WPoint2d wp1 = {{0.0, 1.0}};
  cout << ExactEqual(wp0, wp1) << ", (c.f. 1)" << endl;

  double v = 0; for (int i = 0; i < 10; ++i) { v += 0.1; }
  Point2d p0 = {{ 0.0, v }};
  Point2d p1 = {{0.0, 1.0}};
  cout << ExactEqual(p0, p1) << ", (c.f. 0)" << endl;

  cout << (v == 1.0) << ", (c.f. 0)" << endl;
  cout << bbrcit::almost_equal(v, 1.0) << ", (c.f. 1)" << endl;

  cout << endl;

  return 0; 
}
