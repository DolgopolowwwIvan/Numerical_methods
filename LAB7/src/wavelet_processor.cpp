#include "../include/wavelet_processor.h"
#include "../include/math_constants.h"
#include "../include/signal_operations.h"
#include "../include/signal_transformer.h"
#include <cmath>
#include <stdexcept>

namespace signal_processing
{
    static int wrap_index(int a, int n)
    {
        int r = a % n;
        return (r < 0) ? (r + n) : r;
    }

    wavelet_processor::wavelet_processor(int data_size, wavelet_type type)
    {
        int n = data_size;
        lowpass_filter.assign(n, std::complex<double>(0.0, 0.0));
        highpass_filter.assign(n, std::complex<double>(0.0, 0.0));

        switch (type)
        {
        case wavelet_type::shannon:
        {
            if (n % 4 != 0)
                throw std::runtime_error("размер данных должен быть кратен 4");

            double sqrt2 = 1.0 / std::sqrt(2.0);
            lowpass_filter[0] = std::complex<double>(sqrt2, 0.0);
            highpass_filter[0] = std::complex<double>(sqrt2, 0.0);

            for (int i = 1; i < n; i++)
            {
                double denominator = std::sin(constants::pi * i / n);
                if (std::abs(denominator) < 1e-10)
                {
                    lowpass_filter[i] = std::complex<double>(0.0, 0.0);
                    highpass_filter[i] = std::complex<double>(0.0, 0.0);
                    continue;
                }

                double real_part = std::sqrt(2.0) / n * std::cos(constants::pi * i / n) *
                    std::sin(constants::pi * i / 2.0) / denominator;
                double imag_part = -std::sqrt(2.0) / n * std::sin(constants::pi * i / n) *
                    std::sin(constants::pi * i / 2.0) / denominator;

                std::complex<double> value(real_part, imag_part);
                lowpass_filter[i] = value;
                highpass_filter[i] = std::complex<double>(((i % 2 == 0) ? 1 : -1) * value.real(),
                    ((i % 2 == 0) ? 1 : -1) * value.imag());
            }
            break;
        }
        case wavelet_type::haar:
        {
            double scale = 1.0 / std::sqrt(2.0);
            lowpass_filter[0] = std::complex<double>(scale, 0.0);
            lowpass_filter[1] = std::complex<double>(scale, 0.0);
            highpass_filter[0] = std::complex<double>(scale, 0.0);
            highpass_filter[1] = std::complex<double>(-scale, 0.0);
            break;
        }
        case wavelet_type::daubechies6:
        {
            const double filter_coeffs[6] = {
                0.3326705529500826,
                0.8068915093110928,
                0.4598775021184915,
                -0.1350110200102546,
                -0.08544127388202666,
                0.03522629188570953
            };

            for (int i = 0; i < 6; i++)
            {
                lowpass_filter[i] = std::complex<double>(filter_coeffs[i], 0.0);
            }

            for (int k = 0; k < n; k++)
            {
                int idx = wrap_index(1 - k, n);
                double sign = (k % 2 == 0) ? -1.0 : 1.0;
                highpass_filter[k] = std::complex<double>(sign * lowpass_filter[idx].real(), 0.0);
            }
            break;
        }
        default:
            throw std::runtime_error("неизвестный тип вейвлета");
        }
    }

    void wavelet_processor::build_filter_system(int stages)
    {
        signal_operations operations;
        int n = (int)lowpass_filter.size();

        std::vector<std::vector<std::complex<double>>> low_filters(stages);
        std::vector<std::vector<std::complex<double>>> high_filters(stages);

        low_filters[0] = lowpass_filter;
        high_filters[0] = highpass_filter;

        for (int i = 1; i < stages; i++)
        {
            int element_count = n / (int)std::pow(2.0, i);
            low_filters[i].assign(element_count, std::complex<double>(0.0, 0.0));
            high_filters[i].assign(element_count, std::complex<double>(0.0, 0.0));

            for (int n_idx = 0; n_idx < element_count; n_idx++)
            {
                int max_idx = (int)std::pow(2.0, i);
                for (int k = 0; k < max_idx; k++)
                {
                    low_filters[i][n_idx] += low_filters[0][n_idx + k * n / max_idx];
                    high_filters[i][n_idx] += high_filters[0][n_idx + k * n / max_idx];
                }
            }
        }

        signal_transformer transformer;
        std::vector<std::complex<double>> upsampled_low, upsampled_high;

        decomposition_filters.resize(stages);
        reconstruction_filters.resize(stages);

        decomposition_filters[0] = high_filters[0];
        reconstruction_filters[0] = low_filters[0];

        for (int i = 1; i < stages; i++)
        {
            operations.apply_upsampling(i, low_filters[i], upsampled_low);
            operations.apply_upsampling(i, high_filters[i], upsampled_high);

            transformer.compute_convolution(reconstruction_filters[i - 1], upsampled_high, decomposition_filters[i]);
            transformer.compute_convolution(reconstruction_filters[i - 1], upsampled_low, reconstruction_filters[i]);
        }
    }

    void wavelet_processor::generate_basis_functions(int stage,
        std::vector<std::vector<std::complex<double>>>& wavelet_basis,
        std::vector<std::vector<std::complex<double>>>& scaling_basis)
    {
        signal_operations operations;

        int data_size = (int)lowpass_filter.size();
        int basis_elements = data_size / (int)std::pow(2.0, stage);

        if ((int)reconstruction_filters.size() < stage)
        {
            build_filter_system(stage + 1);
        }

        wavelet_basis.resize(basis_elements);
        scaling_basis.resize(basis_elements);

        for (int i = 0; i < basis_elements; i++)
        {
            int shift_amount = (int)std::pow(2.0, stage) * i;

            std::vector<std::complex<double>> shifted_wavelet;
            operations.perform_circular_shift(shift_amount, decomposition_filters[stage - 1], shifted_wavelet);
            wavelet_basis[i] = shifted_wavelet;

            std::vector<std::complex<double>> shifted_scaling;
            operations.perform_circular_shift(shift_amount, reconstruction_filters[stage - 1], shifted_scaling);
            scaling_basis[i] = shifted_scaling;
        }
    }

    void wavelet_processor::perform_decomposition(int stage,
        const std::vector<std::complex<double>>& input_signal,
        std::vector<std::complex<double>>& wavelet_coeffs,
        std::vector<std::complex<double>>& scaling_coeffs)
    {
        signal_operations operations;

        std::vector<std::vector<std::complex<double>>> wavelet_basis, scaling_basis;
        generate_basis_functions(stage, wavelet_basis, scaling_basis);

        int basis_elements = (int)wavelet_basis.size();
        wavelet_coeffs.assign(basis_elements, std::complex<double>(0.0, 0.0));
        scaling_coeffs.assign(basis_elements, std::complex<double>(0.0, 0.0));

        for (int basis_idx = 0; basis_idx < basis_elements; basis_idx++)
        {
            wavelet_coeffs[basis_idx] = operations.compute_dot_product(input_signal, wavelet_basis[basis_idx]);
            scaling_coeffs[basis_idx] = operations.compute_dot_product(input_signal, scaling_basis[basis_idx]);
        }
    }

    void wavelet_processor::perform_reconstruction(int stage,
        const std::vector<std::complex<double>>& wavelet_coeffs,
        const std::vector<std::complex<double>>& scaling_coeffs,
        std::vector<std::complex<double>>& lowpass_part,
        std::vector<std::complex<double>>& highpass_part,
        std::vector<std::complex<double>>& reconstructed_signal)
    {
        std::vector<std::vector<std::complex<double>>> wavelet_basis, scaling_basis;
        generate_basis_functions(stage, wavelet_basis, scaling_basis);

        int basis_elements = (int)wavelet_basis.size();
        int data_size = (int)lowpass_filter.size();

        lowpass_part.assign(data_size, std::complex<double>(0.0, 0.0));
        highpass_part.assign(data_size, std::complex<double>(0.0, 0.0));
        reconstructed_signal.assign(data_size, std::complex<double>(0.0, 0.0));

        for (int data_idx = 0; data_idx < data_size; data_idx++)
        {
            std::complex<double> lowpass_component(0.0, 0.0);
            std::complex<double> highpass_component(0.0, 0.0);

            for (int basis_idx = 0; basis_idx < basis_elements; basis_idx++)
            {
                lowpass_component = lowpass_component + (scaling_coeffs[basis_idx] * scaling_basis[basis_idx][data_idx]);
                highpass_component = highpass_component + (wavelet_coeffs[basis_idx] * wavelet_basis[basis_idx][data_idx]);
            }

            lowpass_part[data_idx] = lowpass_component;
            highpass_part[data_idx] = highpass_component;
            reconstructed_signal[data_idx] = lowpass_component + highpass_component;
        }
    }
}