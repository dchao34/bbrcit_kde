#ifndef BBRCITKDE_FFT_H__
#define BBRCITKDE_FFT_H__

#include <vector>
#include <complex>

void fft(const std::vector<std::complex<double>> &a, std::vector<std::complex<double>> &y, unsigned r);
void ifft(const std::vector<std::complex<double>> &a, std::vector<std::complex<double>> &y, unsigned r);

#endif
