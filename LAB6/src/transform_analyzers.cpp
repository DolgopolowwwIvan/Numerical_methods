#include "transform_analyzers.h"
#include "complex_operations.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace std::chrono;

AnalysisResults analyzeSignal(const vector<Complex>& signal) {
    AnalysisResults results;

    auto start_dft = high_resolution_clock::now();
    results.dft_result = dft(signal);
    auto end_dft = high_resolution_clock::now();

    auto start_fft = high_resolution_clock::now();
    results.fft_result = fft(signal);
    auto end_fft = high_resolution_clock::now();

    results.timing.dft_time = duration_cast<microseconds>(end_dft - start_dft).count();
    results.timing.fft_time = duration_cast<microseconds>(end_fft - start_fft).count();

    return results;
}

void printResultsTable(const vector<Complex>& signal, const vector<Complex>& dft_result) {
    cout << fixed << setprecision(6);
    cout << setw(4) << "m" << setw(12) << "Re z" << setw(15) << "Re z_hat"
        << setw(15) << "Im z_hat" << setw(15) << "Amplitude" << setw(12) << "Phase" << endl;
    cout << string(75, '-') << endl;

    int N = signal.size();
    double amplitude_limit = 1e-6;
    int count = 0;

    for (int m = 0; m < N; m++) {
        double amplitude = abs(dft_result[m]);
        double phase = arg(dft_result[m]);

        bool significant_amplitude = (amplitude > amplitude_limit);

        if (significant_amplitude) { 
            cout << setw(4) << m << setw(12) << signal[m].real()
                << setw(15) << dft_result[m].real()
                << setw(15) << dft_result[m].imag()
                << setw(15) << amplitude
                << setw(12) << phase << endl;
            count++;
        }
    }

    cout << defaultfloat;
    cout << "Total components output: " << count << " from " << N << endl;
}