#ifndef DATA_EXPORTERS_H
#define DATA_EXPORTERS_H

#include <vector>
#include <complex>
#include <string>

using Complex = std::complex<double>;

// Экспорт данных
void exportToCSV(const std::string& filename,
                 const std::vector<Complex>& original_signal,
                 const std::vector<Complex>& filtered_signal,
                 const std::vector<Complex>& dft_original,
                 const std::vector<Complex>& dft_filtered);

void exportDiscontinuousSignal(const std::string& filename, 
                               const std::vector<Complex>& signal);

#endif