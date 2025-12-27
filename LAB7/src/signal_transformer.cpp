#include "../include/signal_transformer.h"
#include "../include/math_constants.h"
#include <cmath>

namespace signal_processing
{
    void signal_transformer::fast_fourier_transform(const std::vector<std::complex<double>>& input,
        std::vector<std::complex<double>>& output)
    {
        int size = (int)input.size();
        int half_size = size / 2;
        output.assign(size, std::complex<double>(0.0, 0.0));

        std::complex<double> exponent, u_part, v_part;

        for (int m = 0; m < half_size; m++)
        {
            u_part = { 0.0, 0.0 };
            v_part = { 0.0, 0.0 };

            for (int n = 0; n < half_size; n++)
            {
                exponent = { std::cos(-constants::two_pi * m * n / half_size),
                           std::sin(-constants::two_pi * m * n / half_size) };
                u_part += input[2 * n] * exponent;
                v_part += input[2 * n + 1] * exponent;
            }

            exponent = { std::cos(-constants::two_pi * m / size),
                       std::sin(-constants::two_pi * m / size) };
            output[m] = u_part + exponent * v_part;
            output[m + half_size] = u_part - exponent * v_part;
        }
    }

    void signal_transformer::inverse_fast_fourier_transform(const std::vector<std::complex<double>>& input,
        std::vector<std::complex<double>>& output)
    {
        int size = (int)input.size();
        fast_fourier_transform(input, output);

        std::complex<double> temp_value;
        for (int i = 1; i <= size / 2; i++)
        {
            temp_value = output[i];
            output[i] = output[size - i] / double(size);
            output[size - i] = temp_value / double(size);
        }
        output[0] /= double(size);
    }

    void signal_transformer::compute_convolution(const std::vector<std::complex<double>>& vector1,
        const std::vector<std::complex<double>>& vector2,
        std::vector<std::complex<double>>& result)
    {
        int size = (int)vector1.size();
        std::vector<std::complex<double>> intermediate(size);

        result.clear();
        result.resize(size);

        fast_fourier_transform(vector1, result);
        fast_fourier_transform(vector2, intermediate);

        for (int i = 0; i < size; i++)
            intermediate[i] *= result[i];

        inverse_fast_fourier_transform(intermediate, result);
    }
}