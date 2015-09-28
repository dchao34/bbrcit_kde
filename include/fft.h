#ifndef __FFT_H__
#define __FFT_H__

#include <vector>
#include <complex>

void fft(const std::vector<std::complex<double>> &a, std::vector<std::complex<double>> &y, unsigned r);
void ifft(const std::vector<std::complex<double>> &a, std::vector<std::complex<double>> &y, unsigned r);

#endif
