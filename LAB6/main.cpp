#include <iostream>
#include <locale>
#include "signal_generator.h"
#include "transform_analyzers.h"
#include "complex_operations.h"
#include "signal_filters.h"
#include "data_exporters.h"
#include "utilities.h"

using namespace std;

int main() {
    setlocale(LC_ALL, "en_US.UTF-8");
    
    SignalParams params;
    params.N = 1024; // 2 ^ 10
    params.A = 2.26;
    params.B = 0.24;
    params.omega1 = 2;
    params.omega2 = 199;
    params.phi = PI / 6;
    
    printSectionHeader("SECTION 2: SIGNAL GENERATION");
    vector<Complex> signal = generateSignal1(params);
    cout << "Signal generated. N = " << params.N << " points" << endl;
    cout << "Parameters: A = " << params.A 
         << ", B = " << params.B 
         << ", ω1 = " << params.omega1 
         << ", ω2 = " << params.omega2 
         << ", φ = " << params.phi << " rad" << endl;
    
    printSectionHeader("SECTION 3: DFT AND FFT");
    AnalysisResults analysis = analyzeSignal(signal);
    
    cout << "DFT execution time: " << analysis.timing.dft_time << " μs" << endl;
    cout << "FFT execution time: " << analysis.timing.fft_time << " μs" << endl;
    cout << "Speedup factor (DFT/FFT): " 
         << analysis.timing.dft_time / analysis.timing.fft_time << endl;
    
    cout << "\nSignificant spectral components (DFT):" << endl;
    printResultsTable(signal, analysis.dft_result);
    
    printSectionHeader("SECTION 4: NOISE COMPONENT FILTERING");
    vector<Complex> filtered_dft = filterHighFrequencies(analysis.dft_result);
    cout << "High-frequency components filtered." << endl;
    cout << "Original DFT size: " << analysis.dft_result.size() << " components" << endl;
    cout << "Filtered DFT size: " << filtered_dft.size() << " components" << endl;
    

    printSectionHeader("SECTION 5: DATA EXPORT FOR VISUALIZATION");
    vector<Complex> reconstructed = idft(filtered_dft);
    
    exportToCSV("fourier_analysis_results.csv", signal, reconstructed,
                analysis.dft_result, filtered_dft);
    
    cout << "Data exported to 'signal_analysis.csv'" << endl;
    cout << "Files created:" << endl;
    cout << "  - Original signal: " << params.N << " points" << endl;
    cout << "  - Reconstructed signal: " << reconstructed.size() << " points" << endl;
    cout << "  - DFT spectrum: " << analysis.dft_result.size() << " components" << endl;
    cout << "  - Filtered spectrum: " << filtered_dft.size() << " components" << endl;
    
    printSectionHeader("SECTION 6: DISCONTINUOUS SIGNAL ANALYSIS");
    analyzeDiscontinuousSignal(params);
    
    cout << "\n" << string(60, '=') << endl;
    cout << "PROGRAM COMPLETED SUCCESSFULLY" << endl;
    cout << string(60, '=') << endl;
    
    return 0;
}