#include <memory>

#include <Kernel.h>
#include <Kde1d.h>

using namespace std;

Kde1d::Kde1d() : bandwidth_(1.0) {
  pkern_ = shared_ptr<Kernel1d>(new GaussKernel1d());
}

Kde1d::Kde1d(const vector<Point1d> &data, double bw, const Kernel1d &kernel) : 
  bandwidth_(bw), data_(data) {
  pkern_ = kernel.clone();
}

double Kde1d::eval(const Point1d &target) const {

  double val = 0.0;
  for (auto p : data_) {
    auto x = (target - p) / bandwidth_; 
    val += (*pkern_)(x);
  }
  val /= bandwidth_ * data_.size();

  return val;
}
