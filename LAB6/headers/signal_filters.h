#ifndef SIGNAL_FILTERS_H
#define SIGNAL_FILTERS_H

#include <vector>
#include <complex>

using Complex = std::complex<double>;

// Фильтры сигналов
std::vector<Complex> filterHighFrequencies(const std::vector<Complex>& dft_result);

#endif