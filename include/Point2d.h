#ifndef BBRCITKDE_POINT2D_H__
#define BBRCITKDE_POINT2D_H__

#ifdef __CUDACC__
#define CUDA_ALIGN16 __align__(16)
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_ALIGN16
#define CUDA_CALLABLE
#endif

namespace bbrcit {

template <typename T>
class Point2d {
  using FloatType = float;
};

template <>
class CUDA_ALIGN16 Point2d<float> {

  public:
    using FloatType = float;

    CUDA_CALLABLE Point2d() = default;
    CUDA_CALLABLE Point2d(float x, float y, float m) : 
      x_(x), y_(y), m_(m), abw_(1.0f) {}
    CUDA_CALLABLE Point2d(float x, float y, float m, float abw) : 
      x_(x), y_(y), m_(m), abw_(abw) {}

    template<typename HostDecPointT>
      Point2d<float>& operator=(const HostDecPointT &pt) {
        x_ = pt[0]; 
        y_ = pt[1];
        m_ = pt.attributes().mass();
        abw_ = pt.attributes().abw();
        return *this;
      }

    CUDA_CALLABLE ~Point2d() = default;
    CUDA_CALLABLE Point2d(const Point2d<float>&) = default;
    CUDA_CALLABLE Point2d<float>& operator=(const Point2d<float>&) = default;

    CUDA_CALLABLE float x() const { return x_; }
    CUDA_CALLABLE float y() const { return y_; }
    CUDA_CALLABLE float m() const { return m_; }
    CUDA_CALLABLE float abw() const { return abw_; }

  private:
    float x_ = 0.0f;
    float y_ = 0.0f;
    float m_ = 1.0f;
    float abw_ = 1.0f;
};

template <>
class CUDA_ALIGN16 Point2d<double> {

  public:
    using FloatType = double;

    CUDA_CALLABLE Point2d() = default;
    CUDA_CALLABLE Point2d(double x, double y, double m) : 
      x_(x), y_(y), m_(m), abw_(1.0) {}
    CUDA_CALLABLE Point2d(double x, double y, double m, double abw) : 
      x_(x), y_(y), m_(m), abw_(abw) {}

    template<typename HostDecPointT>
      Point2d<double>& operator=(const HostDecPointT &pt) {
        x_ = pt[0]; 
        y_ = pt[1];
        m_ = pt.attributes().mass();
        abw_ = pt.attributes().abw();
        return *this;
      }

    CUDA_CALLABLE ~Point2d() = default;
    CUDA_CALLABLE Point2d(const Point2d<double>&) = default;
    CUDA_CALLABLE Point2d<double>& operator=(const Point2d<double>&) = default;

    CUDA_CALLABLE double x() const { return x_; }
    CUDA_CALLABLE double y() const { return y_; }
    CUDA_CALLABLE double m() const { return m_; }
    CUDA_CALLABLE double abw() const { return abw_; }

  private:
    double x_ = 0.0;
    double y_ = 0.0;
    double m_ = 1.0;
    double abw_ = 1.0;
};


}

#endif
