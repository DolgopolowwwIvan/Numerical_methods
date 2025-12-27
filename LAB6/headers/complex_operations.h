#ifndef COMPLEX_OPERATIONS_H
#define COMPLEX_OPERATIONS_H

#include <vector>
#include <complex>

using Complex = std::complex<double>;
const double PI = 3.14159265358979323846;

std::vector<Complex> dft(const std::vector<Complex>& input);
std::vector<Complex> idft(const std::vector<Complex>& input);
std::vector<Complex> fft(const std::vector<Complex>& input);
std::vector<Complex> ifft(const std::vector<Complex>& input);

extern const double PI;

#endif