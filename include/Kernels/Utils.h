#ifndef __KERNEL_UTILS_H__
#define __KERNEL_UTILS_H__

template <typename T>
inline T epanechnikov_choice(T v1, T v2, T v3) {
  return (std::abs(v3) >= std::abs(v2) && std::abs(v3) >= std::abs(v1)) ? v2 : v3;
}

#endif
