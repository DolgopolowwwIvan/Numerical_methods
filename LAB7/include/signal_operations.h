#pragma once
#ifndef SIGNAL_OPERATIONS_H
#define SIGNAL_OPERATIONS_H

#include <complex>
#include <vector>

namespace signal_processing
{
    class signal_operations
    {
    public:
        void perform_circular_shift(int shift_amount,
            const std::vector<std::complex<double>>& data,
            std::vector<std::complex<double>>& result);

        void apply_downsampling(int level,
            const std::vector<std::complex<double>>& data,
            std::vector<std::complex<double>>& result);

        void apply_upsampling(int level,
            const std::vector<std::complex<double>>& data,
            std::vector<std::complex<double>>& result);

        std::complex<double> compute_dot_product(const std::vector<std::complex<double>>& vec1,
            const std::vector<std::complex<double>>& vec2);
    };
}

#endif