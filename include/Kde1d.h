#include <Point.h>
#include <Kernel.h>

#include <vector>
#include <memory>

class Kde1d {

  public:
    Kde1d();
    ~Kde1d() = default;
    Kde1d(const std::vector<Point1d> &data, double bw=1.0, const Kernel1d &kernel=GaussKernel1d());

    double get_bandwidth() const { return bandwidth_; }

    void set_bandwidth(double bw) { bandwidth_ = bw; }
    void set_kernel(const Kernel1d &kernel) { pkern_ = kernel.clone(); }
    void set_data(const std::vector<Point1d> &data) { data_ = data; }

    double eval(const Point1d&) const;

  private:
    double bandwidth_;
    std::shared_ptr<Kernel1d> pkern_;
    std::vector<Point1d> data_;

};
