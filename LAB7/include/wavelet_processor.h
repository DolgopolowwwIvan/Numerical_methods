#pragma once
#ifndef WAVELET_PROCESSOR_H
#define WAVELET_PROCESSOR_H

#include <vector>
#include <complex>
#include "signal_operations.h"
#include "signal_transformer.h"

namespace signal_processing
{
    class wavelet_processor
    {
    private:
        std::vector<std::complex<double>> lowpass_filter, highpass_filter;
        std::vector<std::vector<std::complex<double>>> decomposition_filters, reconstruction_filters;

    public:
        enum class wavelet_type
        {
            haar = 1,
            shannon = 2,
            daubechies6 = 3
        };

        wavelet_processor(int data_size, wavelet_type type);

    private:
        void build_filter_system(int stages);

        void generate_basis_functions(int stage,
            std::vector<std::vector<std::complex<double>>>& wavelet_basis,
            std::vector<std::vector<std::complex<double>>>& scaling_basis);

    public:
        void perform_decomposition(int stage,
            const std::vector<std::complex<double>>& input_signal,
            std::vector<std::complex<double>>& wavelet_coeffs,
            std::vector<std::complex<double>>& scaling_coeffs);

        void perform_reconstruction(int stage,
            const std::vector<std::complex<double>>& wavelet_coeffs,
            const std::vector<std::complex<double>>& scaling_coeffs,
            std::vector<std::complex<double>>& lowpass_part,
            std::vector<std::complex<double>>& highpass_part,
            std::vector<std::complex<double>>& reconstructed_signal);
    };
}

#endif