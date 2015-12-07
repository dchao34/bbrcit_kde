#include <stdexcept>
#include <Point.h>

Point1d::Point1d(double s) : x_(s) {}

Point1d& Point1d::operator*=(double c) { x_ *= c; return *this; }
Point1d& Point1d::operator/=(double c) { x_ /= c; return *this; }
Point1d& Point1d::operator+=(const Point1d &p) { x_ += p.x_; return *this; }
Point1d& Point1d::operator-=(const Point1d &p) { x_ -= p.x_; return *this; }

const double& Point1d::operator[](std::size_t i) const {
  if (i != 0) { throw std::out_of_range("Point1d: operator[] accessed with non-zero index."); }
  return x_;
}

double& Point1d::operator[](std::size_t i) {
  return const_cast<double&>(static_cast<const Point1d&>(*this)[i]);
}

Point1d operator+(const Point1d &lhs, const Point1d &rhs) {
  Point1d result(lhs); result += rhs; return result;
}

Point1d operator-(const Point1d &lhs, const Point1d &rhs) {
  Point1d result(lhs); result -= rhs; return result;
}

Point1d operator*(double c, const Point1d &p) { return p * c; }
Point1d operator*(const Point1d &p, double c) {
  Point1d result(p); result *= c; return result;
}

Point1d operator/(const Point1d &p, double c) {
  Point1d result(p); result /= c; return result;
}

bool operator==(const Point1d &lhs, const Point1d &rhs) { return lhs.x() == rhs.x(); }
bool operator!=(const Point1d &lhs, const Point1d &rhs) { return !(lhs == rhs); }
bool operator<(const Point1d &lhs, const Point1d &rhs) { return lhs.x() < rhs.x(); }
bool operator<=(const Point1d &lhs, const Point1d &rhs) { return !(rhs < lhs); }
bool operator>(const Point1d &lhs, const Point1d &rhs) { return rhs < lhs; }
bool operator>=(const Point1d &lhs, const Point1d &rhs) { return !(lhs < rhs); }
