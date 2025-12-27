#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "../include/wavelet_processor.h"
#include "../include/math_constants.h"

static void save_signal_to_csv(const std::string& filepath, const std::vector<std::complex<double>>& signal)
{
    std::ofstream file(filepath);
    file << "index,real_part,imag_part\n";
    for (int i = 0; i < (int)signal.size(); i++)
        file << i << "," << std::setprecision(17) << signal[i].real() << "," << signal[i].imag() << "\n";
}

static void save_coefficients_to_csv(const std::string& filepath, int stage,
    const std::vector<std::complex<double>>& psi_coeffs,
    const std::vector<std::complex<double>>& phi_coeffs)
{
    std::ofstream file(filepath);
    file << "k,index,psi_real,psi_imag,phi_real,phi_imag,psi_magnitude,phi_magnitude\n";
    for (int k = 0; k < (int)psi_coeffs.size(); k++)
    {
        int idx = (int)(std::pow(2.0, stage) * k);
        file << k << "," << idx << ","
            << std::setprecision(17) << psi_coeffs[k].real() << "," << psi_coeffs[k].imag() << ","
            << phi_coeffs[k].real() << "," << phi_coeffs[k].imag() << ","
            << std::abs(psi_coeffs[k]) << "," << std::abs(phi_coeffs[k]) << "\n";
    }
}

static void save_filter_results_to_csv(const std::string& filepath,
    const std::vector<std::complex<double>>& original,
    const std::vector<std::complex<double>>& filtered,
    const std::vector<std::complex<double>>& difference)
{
    std::ofstream file(filepath);
    file << "index,original_real,filtered_real,difference_real\n";
    for (int i = 0; i < (int)original.size(); i++)
        file << i << "," << std::setprecision(17)
        << original[i].real() << ","
        << filtered[i].real() << ","
        << difference[i].real() << "\n";
}

static void save_pq_components_to_csv(const std::string& filepath,
    const std::vector<std::complex<double>>& p_component,
    const std::vector<std::complex<double>>& q_component)
{
    std::ofstream file(filepath);
    file << "index,P_real,Q_real\n";
    for (int i = 0; i < (int)p_component.size(); i++)
        file << i << "," << std::setprecision(17)
        << p_component[i].real() << ","
        << q_component[i].real() << "\n";
}

static std::string get_wavelet_name(signal_processing::wavelet_processor::wavelet_type type)
{
    using wt = signal_processing::wavelet_processor::wavelet_type;
    if (type == wt::haar) return "haar";
    if (type == wt::shannon) return "shannon";
    return "d6";
}

static void compute_only_p_component(signal_processing::wavelet_processor& processor,
    int stage,
    const std::vector<std::complex<double>>& signal,
    std::vector<std::complex<double>>& p_component)
{
    std::vector<std::complex<double>> psi_coeffs, phi_coeffs;
    processor.perform_decomposition(stage, signal, psi_coeffs, phi_coeffs);

    std::vector<std::complex<double>> zeroed_psi(psi_coeffs.size(), { 0.0, 0.0 });
    std::vector<std::complex<double>> temp_p, temp_q, temp_recovery;
    processor.perform_reconstruction(stage, zeroed_psi, phi_coeffs, temp_p, temp_q, temp_recovery);
    p_component = temp_p;
}

static void process_wavelet_basis(const std::string& output_directory,
    signal_processing::wavelet_processor::wavelet_type type,
    const std::vector<std::complex<double>>& input_signal,
    int max_stages)
{
    int n = (int)input_signal.size();
    signal_processing::wavelet_processor processor(n, type);
    std::string basis_name = get_wavelet_name(type);

    // Создаем директорию, если она не существует
    std::string cmd = "mkdir -p " + output_directory;
    system(cmd.c_str());

    for (int level = 1; level <= max_stages; level++)
    {
        std::vector<std::complex<double>> psi_coeffs, phi_coeffs;
        processor.perform_decomposition(level, input_signal, psi_coeffs, phi_coeffs);

        save_coefficients_to_csv(output_directory + "/coeffs_before_" + basis_name + "_stage" + std::to_string(level) + ".csv",
            level, psi_coeffs, phi_coeffs);

        std::vector<std::complex<double>> zeroed_psi(psi_coeffs.size(), { 0.0, 0.0 });
        save_coefficients_to_csv(output_directory + "/coeffs_after_" + basis_name + "_stage" + std::to_string(level) + ".csv",
            level, zeroed_psi, phi_coeffs);

        std::vector<std::complex<double>> p_comp, q_comp, filtered_signal;
        processor.perform_reconstruction(level, zeroed_psi, phi_coeffs, p_comp, q_comp, filtered_signal);

        std::vector<std::complex<double>> difference(n);
        for (int i = 0; i < n; i++)
            difference[i] = input_signal[i] - filtered_signal[i];

        save_filter_results_to_csv(output_directory + "/filter_results_" + basis_name + "_stage" + std::to_string(level) + ".csv",
            input_signal, filtered_signal, difference);

        std::vector<std::complex<double>> previous_p(n), previous_q(n);
        if (level == 1)
        {
            previous_p = filtered_signal;
            for (int i = 0; i < n; i++)
                previous_q[i] = { 0.0, 0.0 };
        }
        else
        {
            compute_only_p_component(processor, level - 1, filtered_signal, previous_p);
            for (int i = 0; i < n; i++)
                previous_q[i] = filtered_signal[i] - previous_p[i];
        }

        save_pq_components_to_csv(output_directory + "/pq_components_" + basis_name + "_stage" + std::to_string(level) + ".csv",
            previous_p, previous_q);
    }

}

// Функция для генерации кусочно-постоянного сигнала (ваш вариант 1)
std::vector<std::complex<double>> generate_piecewise_constant_signal(size_t n, double a, double b, double w2)
{
    std::vector<std::complex<double>> signal(n, { 0.0, 0.0 });

    //Используем современный ГСЧ для добавления шума
    // unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    // std::mt19937 gen(seed);
    // std::uniform_real_distribution<double> noise_dist(-0.05, 0.05);

    const double quarter = n / 4.0;
    const double half = n / 2.0;
    const double three_quarters = 3.0 * n / 4.0;

    for (size_t j = 0; j < n; ++j) {
        double value = 0.0;
        
        if ((j >= quarter && j <= half) || (j > three_quarters)) {
            // На участках где сигнал ненулевой: A + B*cos(2πω₂j/N)
            value = a + b * std::cos(2.0 * signal_processing::constants::pi * w2 * j / n);
        }
        // На остальных участках value остаётся 0.0
        
        // Добавляем небольшой шум
        //double noise = noise_dist(gen);
        signal[j] = { value, 0.0 };
    }

    return signal;
}

int main() {
    using namespace signal_processing;
    const int n_power = 9;
    const int n_size = 1 << n_power;  // 2^9 = 512

    const double a_val = 2.87;      // Ваш A = 2.87
    const double b_val = 0.26;      // Ваш B = 0.26
    const int w1_val = 2;           // Ваш ω₁ = 2 (не используется в кусочно-постоянном сигнале)
    const int w2_val = 194;         // Ваш ω₂ = 194
    const int max_stages = 4;

    const std::string output_directory = "./output";

    // Генерация кусочно-постоянного сигнала (вариант 1 из задания)
    std::vector<std::complex<double>> signal_data = 
        generate_piecewise_constant_signal(n_size, a_val, b_val, w2_val);

    save_signal_to_csv(output_directory + "/generated_signal.csv", signal_data);

    // Обработка для всех трех базисов
    process_wavelet_basis(output_directory,
        wavelet_processor::wavelet_type::haar,
        signal_data, max_stages);

    process_wavelet_basis(output_directory,
        wavelet_processor::wavelet_type::shannon,
        signal_data, max_stages);

    process_wavelet_basis(output_directory,
        wavelet_processor::wavelet_type::daubechies6,
        signal_data, max_stages);

    std::cout << "files saved in directory: " << output_directory << std::endl;

    return 0;
}