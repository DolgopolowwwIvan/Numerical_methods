#include "signal_filters.h"
#include <iostream>
#include <cmath>

using namespace std;

vector<Complex> filterHighFrequencies(const vector<Complex>& dft_result) {
    int N = dft_result.size();
    vector<Complex> filtered = dft_result;

    int keep_count = N / 10; 

    cout << "Simple filtering: " << endl;
    cout << "We save frequencies: m = 0..." << keep_count << " and " << N - keep_count << "...N-1" << endl;

    int removed_count = 0;
    for (int k = keep_count + 1; k < N - keep_count; k++) {
        if (abs(filtered[k]) > 1e-6) {
            removed_count++;
        }
        filtered[k] = 0;
    }

    cout << "Removed component with non-zero amplitude: " << removed_count << endl;
    return filtered;
}