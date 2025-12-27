#include "../include/signal_operations.h"
#include <cmath>

namespace signal_processing
{
    void signal_operations::perform_circular_shift(int shift_amount,
        const std::vector<std::complex<double>>& data,
        std::vector<std::complex<double>>& result)
    {
        int size = (int)data.size();
        result.assign(size, std::complex<double>(0.0, 0.0));

        for (int i = 0; i < size; i++)
        {
            int shifted_index = i - shift_amount;
            if (shifted_index < 0) shifted_index += size;
            result[i] = data[shifted_index];
        }
    }

    void signal_operations::apply_downsampling(int level,
        const std::vector<std::complex<double>>& data,
        std::vector<std::complex<double>>& result)
    {
        int factor = (int)std::pow(2.0, level);
        int new_size = (int)data.size() / factor;
        result.assign(new_size, std::complex<double>(0.0, 0.0));

        for (int i = 0; i < new_size; i++)
            result[i] = data[i * factor];
    }

    void signal_operations::apply_upsampling(int level,
        const std::vector<std::complex<double>>& data,
        std::vector<std::complex<double>>& result)
    {
        int factor = (int)std::pow(2.0, level);
        int new_size = (int)data.size() * factor;
        result.assign(new_size, std::complex<double>(0.0, 0.0));

        for (int i = 0; i < new_size; i++)
        {
            if (i % factor == 0)
                result[i] = data[i / factor];
            else
                result[i] = std::complex<double>(0.0, 0.0);
        }
    }

    std::complex<double> signal_operations::compute_dot_product(const std::vector<std::complex<double>>& vec1,
        const std::vector<std::complex<double>>& vec2)
    {
        int size = (int)vec1.size();
        std::complex<double> result(0.0, 0.0);

        for (int i = 0; i < size; i++)
            result += vec1[i] * std::conj(vec2[i]);

        return result;
    }
}