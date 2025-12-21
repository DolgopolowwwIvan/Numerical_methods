#ifndef TRANSFORM_ANALYZERS_H
#define TRANSFORM_ANALYZERS_H

#include <vector>
#include <complex>
#include <chrono>

using Complex = std::complex<double>;

struct TimingResults{
    long long dft_time;
    long long fft_time;
};

struct AnalysisResults{
    std::vector<Complex> dft_result;
    std::vector<Complex> fft_result;
    TimingResults timing;
};

AnalysisResults analyzeSignal(const std::vector<Complex>& signal);
void printResultsTable(const std::vector<Complex>& signal, 
                       const std::vector<Complex>& dft_result);

#endif