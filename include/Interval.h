#ifndef BBRCITKDE_INTERVAL_H__
#define BBRCITKDE_INTERVAL_H__

#include <iostream>
#include <cmath>
#include <exception>

// API
// ---

namespace bbrcit {

template <typename T> class Interval;

template <typename T> 
std::ostream& operator<<(std::ostream&, const Interval<T>&);

// test whether two Interval<>'s intersect.
template <typename T> 
bool intersect(const Interval<T>&, const Interval<T>&);

template <typename T=double>
class Interval {

  public:
    Interval();
    Interval(const T&, const T&);
    Interval(const Interval<T>&) = default;
    Interval(Interval<T>&&) = default;
    Interval<T>& operator=(const Interval<T>&) = default;
    Interval<T>& operator=(Interval<T>&&) = default;
    ~Interval();

    // returns the length
    T length() const;

    // returns the lower/upper endpoints and the midpoint
    T lower() const;
    T middle() const;
    T upper() const;

    // set new endpoint values for this Interval<>. 
    void resize(const T&, const T&);

    // returns true if the argument is fully contained in this Interval<>.
    bool contains(const T&) const;
    bool contains(const Interval<T>&) const;

    // returns min/max distance argument to any point in this Interval<>.
    T min_dist(const T&) const;
    T max_dist(const T&) const;

  private:
    T lower_, upper_;
};


// Implementations
// ---------------

template <typename T> 
std::ostream& operator<<(std::ostream &os, const Interval<T> &i) {
  os << "(" << i.lower() << ", " << i.upper() << ")";
  return os;
}

template <typename T>
Interval<T>::Interval() : lower_(T()), upper_(T()) {}

template <typename T>
Interval<T>::Interval(const T &l, const T &u) : lower_(l), upper_(u) {
  using std::swap;
  if (upper_ < lower_) { swap(lower_, upper_); }
}

template <typename T>
Interval<T>::~Interval() {}

template <typename T>
inline T Interval<T>::length() const { return upper_ - lower_; }

template <typename T>
inline T Interval<T>::middle() const { return lower_ + (upper_ - lower_) / 2; }

template <typename T>
inline T Interval<T>::lower() const { return lower_; }

template <typename T>
inline T Interval<T>::upper() const { return upper_; }

template <typename T> 
inline bool intersect(const Interval<T> &i1, const Interval<T> &i2) {
  return i1.lower() <= i2.upper() && i2.lower() <= i1.upper();
}

template <typename T> 
inline bool Interval<T>::contains(const T &p) const {
  return lower_ <= p && p <= upper_;
}

template <typename T> 
inline bool Interval<T>::contains(const Interval<T> &rhs) const {
  return lower_ <= rhs.lower_ && rhs.upper_ <= upper_;
}

template <typename T> 
T Interval<T>::min_dist(const T &p) const {
  if (this->contains(p)) { return 0; }
  T lower_dist = std::abs(p - lower_);
  T upper_dist = std::abs(p - upper_);
  return lower_dist < upper_dist ? lower_dist : upper_dist;
}

template <typename T> 
T Interval<T>::max_dist(const T &p) const {
  T lower_dist = std::abs(p - lower_);
  T upper_dist = std::abs(p - upper_);
  return lower_dist < upper_dist ? upper_dist : lower_dist;
}

}


#endif
