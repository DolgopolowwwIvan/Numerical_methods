#include "signal_generator.h"
#include "complex_operations.h"
#include <cmath>

using namespace std;

vector<Complex> generateSignal1(const SignalParams& params) {
    vector<Complex> signal(params.N);
    for (int j = 0; j < params.N; j++) {
        double value = params.A * cos(2 * PI * params.omega1 * j / params.N + params.phi) +
            params.B * cos(2 * PI * params.omega2 * j / params.N);
        signal[j] = value;
    }
    return signal;
}

vector<Complex> generateSignal2(const SignalParams& params) {
    vector<Complex> signal(params.N, 0);
    for (int j = params.N / 4; j <= params.N / 2; j++) {
        signal[j] = params.A + params.B * cos(2 * PI * params.omega2 * j / params.N);
    }
    for (int j = 3 * params.N / 4; j < params.N; j++) {
        signal[j] = params.A + params.B * cos(2 * PI * params.omega2 * j / params.N);
    }
    return signal;
}