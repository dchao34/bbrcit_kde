#include <cmath>

class Point {
  public:

    Point() {};
    virtual ~Point() {};

    virtual double norm2() const = 0;
};

class Point1d : public Point {
  public:

    Point1d() = default;
    Point1d(double);

    ~Point1d() = default;
    Point1d(const Point1d&) = default;
    Point1d& operator=(const Point1d&) = default;

    // scaler multiplication and division
    Point1d& operator*=(double);
    Point1d& operator/=(double);

    // vector addition and subtraction
    Point1d& operator+=(const Point1d&);
    Point1d& operator-=(const Point1d&);

    // subscript: index 0 returns reference to x_.
    const double& operator[](std::size_t) const;
    double& operator[](std::size_t);

    // access coordinate value
    double x() const { return x_; }

    // L2 norm and its square
    double norm() const { return std::abs(x_); }
    double norm2() const { return x_*x_; }

  private:
    double x_ = 0.0;
};

Point1d operator*(const Point1d&, double);
Point1d operator*(double, const Point1d&);
Point1d operator/(const Point1d&, double);
Point1d operator+(const Point1d&, const Point1d&);
Point1d operator-(const Point1d&, const Point1d&);
bool operator==(const Point1d&, const Point1d&);
bool operator!=(const Point1d&, const Point1d&);
bool operator<(const Point1d&, const Point1d&);
bool operator<=(const Point1d&, const Point1d&);
bool operator>(const Point1d&, const Point1d&);
bool operator>=(const Point1d&, const Point1d&);
