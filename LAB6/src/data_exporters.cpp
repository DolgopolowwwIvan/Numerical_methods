#include "data_exporters.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

void exportToCSV(const string& filename,
                 const vector<Complex>& original_signal,
                 const vector<Complex>& filtered_signal,
                 const vector<Complex>& dft_original,
                 const vector<Complex>& dft_filtered) {
    ofstream file(filename);
    
    file << "sample,original,filtered,"
         << "spectrum_real,spectrum_imag,spectrum_mag,"
         << "spectrum_filt_real,spectrum_filt_imag,spectrum_filt_mag" << endl;
    
    int N = original_signal.size();
    for (int i = 0; i < N; i++) {
        file << i << ","
             << original_signal[i].real() << ","
             << filtered_signal[i].real() << ","
             << dft_original[i].real() << ","
             << dft_original[i].imag() << ","
             << abs(dft_original[i]) << ","
             << dft_filtered[i].real() << ","
             << dft_filtered[i].imag() << ","
             << abs(dft_filtered[i]) << endl;
    }
    
    file.close();
    cout << "Data exported to: " << filename << endl;
}

void exportDiscontinuousSignal(const string& filename, 
                               const vector<Complex>& signal) {
    ofstream file(filename);
    
    file << "sample,discontinuous_signal" << endl;
    
    int N = signal.size();
    for (int i = 0; i < N; i++) {
        file << i << "," << signal[i].real() << endl;
    }
    
    file.close();
    cout << "Discontinuous signal exported to: " << filename << endl;
}