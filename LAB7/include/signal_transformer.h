#pragma once
#ifndef SIGNAL_TRANSFORMER_H
#define SIGNAL_TRANSFORMER_H

#include <vector>
#include <complex>

namespace signal_processing
{
    class signal_transformer
    {
    public:
        void fast_fourier_transform(const std::vector<std::complex<double>>& input,
            std::vector<std::complex<double>>& output);

        void inverse_fast_fourier_transform(const std::vector<std::complex<double>>& input,
            std::vector<std::complex<double>>& output);

        void compute_convolution(const std::vector<std::complex<double>>& vector1,
            const std::vector<std::complex<double>>& vector2,
            std::vector<std::complex<double>>& result);
    };
}

#endif