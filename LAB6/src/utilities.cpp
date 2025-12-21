#include "utilities.h"
#include "signal_generator.h"
#include "data_exporters.h"
#include <iostream>

using namespace std; 

void printSectionHeader(const std::string& title) {
    std::cout << "\n" << std::string(50, '=') << "\n";
    std::cout << title << "\n\n";
}

void analyzeDiscontinuousSignal(const SignalParams& params) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "POIN 6" << "\n";
    
    std::vector<Complex> signal = generateSignal2(params);
    exportDiscontinuousSignal("signal_with_discontinuities.csv", signal);
}