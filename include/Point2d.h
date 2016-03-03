#ifndef BBRCITKDE_POINT2D_H__
#define BBRCITKDE_POINT2D_H__

#ifdef __CUDACC__
#define CUDA_ALIGN16 __align__(16)
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_ALIGN16
#define CUDA_CALLABLE
#endif

template <typename T>
class Point2d {
  using FloatType = float;
};

template <>
class CUDA_ALIGN16 Point2d<float> {

  public:
    using FloatType = float;

    CUDA_CALLABLE Point2d() = default;
    CUDA_CALLABLE Point2d(float, float, float);

    CUDA_CALLABLE ~Point2d() = default;
    CUDA_CALLABLE Point2d(const Point2d<float>&) = default;
    CUDA_CALLABLE Point2d<float>& operator=(const Point2d<float>&) = default;

    CUDA_CALLABLE float x() const { return x_; }
    CUDA_CALLABLE float y() const { return y_; }
    CUDA_CALLABLE float w() const { return w_; }

  private:
    float x_ = 0.0f;
    float y_ = 0.0f;
    float w_ = 0.0f;
};

inline
Point2d<float>::Point2d(float x, float y, float w) : x_(x), y_(y), w_(w) {}

#endif
