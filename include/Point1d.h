#ifndef BBRCITKDE_POINT1D_H__
#define BBRCITKDE_POINT1D_H__

#ifdef __CUDACC__
#define CUDA_ALIGN16 __align__(16)
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_ALIGN16
#define CUDA_CALLABLE
#endif

namespace bbrcit {

template <typename T>
class Point1d {
  using FloatType = float;
};

template <>
class CUDA_ALIGN16 Point1d<float> {

  public:
    using FloatType = float;

    CUDA_CALLABLE Point1d() = default;
    CUDA_CALLABLE Point1d(float x, float m) : 
      x_(x), m_(m), abw_(1.0f) {}
    CUDA_CALLABLE Point1d(float x, float m, float abw) : 
      x_(x), m_(m), abw_(abw) {}

    template<typename HostDecPointT>
      Point1d<float>& operator=(const HostDecPointT &pt) {
        x_ = pt[0]; 
        m_ = pt.attributes().mass();
        abw_ = pt.attributes().abw();
        return *this;
      }

    CUDA_CALLABLE ~Point1d() = default;
    CUDA_CALLABLE Point1d(const Point1d<float>&) = default;
    CUDA_CALLABLE Point1d<float>& operator=(const Point1d<float>&) = default;

    CUDA_CALLABLE float x() const { return x_; }
    CUDA_CALLABLE float m() const { return m_; }
    CUDA_CALLABLE float abw() const { return abw_; }

  private:
    float x_ = 0.0f;
    float m_ = 1.0f;
    float abw_ = 1.0f;
};

template <>
class CUDA_ALIGN16 Point1d<double> {

  public:
    using FloatType = double;

    CUDA_CALLABLE Point1d() = default;
    CUDA_CALLABLE Point1d(double x, double m) : 
      x_(x), m_(m), abw_(1.0f) {}
    CUDA_CALLABLE Point1d(double x, double m, double abw) : 
      x_(x), m_(m), abw_(abw) {}

    template<typename HostDecPointT>
      Point1d<double>& operator=(const HostDecPointT &pt) {
        x_ = pt[0]; 
        m_ = pt.attributes().mass();
        abw_ = pt.attributes().abw();
        return *this;
      }

    CUDA_CALLABLE ~Point1d() = default;
    CUDA_CALLABLE Point1d(const Point1d<double>&) = default;
    CUDA_CALLABLE Point1d<double>& operator=(const Point1d<double>&) = default;

    CUDA_CALLABLE double x() const { return x_; }
    CUDA_CALLABLE double m() const { return m_; }
    CUDA_CALLABLE double abw() const { return abw_; }

  private:
    double x_ = 0.0f;
    double m_ = 1.0f;
    double abw_ = 1.0f;
};


}

#endif
