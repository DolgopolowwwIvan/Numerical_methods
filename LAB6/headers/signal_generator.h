#ifndef SIGNAL_GENERATOR_H
#define SIGNAL_GENERATOR_H

#include <vector>
#include <complex>

using Complex = std::complex<double>;

struct SignalParams {
    int N;
    double A, B;
    double omega1, omega2;
    double phi;
};

std::vector<Complex> generateSignal1(const SignalParams& params);
std::vector<Complex> generateSignal2(const SignalParams& params);

#endif