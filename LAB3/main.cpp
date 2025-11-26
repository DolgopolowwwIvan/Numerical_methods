#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include "Smoothing_Spline_1D.h"


using namespace Com_Methods;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


std::string get_output_path(const std::string& filename) {
    return "E:/cm3/output/" + filename;
}

void generate_random_data(int N, double M, double sigma,
    std::vector<Point>& points,
    std::vector<double>& values) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(M, sigma);

    points.clear();
    values.clear();

    for (int i = 0; i < N; ++i) {
        double x = i;
        double y = dist(gen);

        points.emplace_back(x, 0, 0);
        values.push_back(y);
    }

    std::string filename = get_output_path("random_data.csv");
    std::ofstream file(filename);

    if (!file.is_open()) {
        return;
    }

    for (int i = 0; i < N; ++i) {
        file << i << "," << values[i] << "\n";
    }
    file.close();

}

// Функция для вычисления статистики
void calculate_statistics(const std::vector<double>& values,
    double& mean, double& min_val,
    double& max_val, double& std_dev) {
    double sum = 0.0;
    double sum_sq = 0.0;
    min_val = values[0];
    max_val = values[0];

    for (const auto& val : values) {
        sum += val;
        sum_sq += val * val;
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
    }

    mean = sum / values.size();
    std_dev = std::sqrt(sum_sq / values.size() - mean * mean);
}

void save_spline_comparison(const std::vector<Point>& points,
    const std::vector<double>& values,
    const std::vector<double>& p_values) {
    std::cout << "Saving spline comparison data for " << points.size() << " points..." << std::endl;

    std::string filename = get_output_path("spline_comparison.csv");
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cout << "ERROR: Cannot open " << filename << std::endl;
        return;
    }

    file << "Index,Original";
    for (double p : p_values) {
        file << ",p=" << p;
    }
    file << ",Interpolation\n";

    std::vector<Smoothing_Spline_1D> smoothing_splines;
    for (double p : p_values) {
        std::cout << "Building smoothing spline with p=" << p << "..." << std::endl;
        Smoothing_Spline_1D spline(p);
        spline.Update_Spline(points, values);
        smoothing_splines.push_back(spline);
    }

    int num_points = points.size();
    for (int i = 0; i < num_points; ++i) {
        double x = points[i].x();
        file << x << "," << values[i];

        for (auto& spline : smoothing_splines) {
            Point p(x, 0, 0);
            double result[3];
            spline.Get_Value(p, result);
            file << "," << result[0];
        }

        Point p(x, 0, 0);
        double result[3];
        file << "," << result[0] << "\n";

        if (i % 100 == 0) {
            std::cout << "Processed " << i << " of " << num_points << " points..." << std::endl;
        }
    }

    file.close();
}

void analyze_smoothing_parameter(const std::vector<Point>& points,
    const std::vector<double>& values,
    const std::vector<double>& p_values) {
    std::cout << "\n--- Analysis of smoothing parameter p influence ---" << std::endl;

    std::string filename = get_output_path("smoothing_analysis.csv");
    std::ofstream analysis_file(filename);

    if (!analysis_file.is_open()) {
        std::cout << "ERROR: Cannot open " << filename << std::endl;
        return;
    }

    analysis_file << "p_value,Average_Deviation,Max_Deviation\n";

    for (double p : p_values) {
        std::cout << "Analyzing p=" << p << "..." << std::endl;
        Smoothing_Spline_1D smoothing_spline(p);
        smoothing_spline.Update_Spline(points, values);

        double total_deviation = 0.0;
        double max_deviation = 0.0;
        int count = points.size();

        for (int i = 0; i < count; i++) {
            Point test_point(points[i].x(), 0, 0);
            double result[3];
            smoothing_spline.Get_Value(test_point, result);

            double deviation = std::abs(result[0] - values[i]);
            total_deviation += deviation;
            if (deviation > max_deviation) {
                max_deviation = deviation;
            }
        }

        double avg_deviation = total_deviation / count;
        analysis_file << p << "," << avg_deviation << "," << max_deviation << "\n";

    }

    analysis_file.close();
}

int main() {

    const int N = 1670;
    const double M = 1.04;
    const double sigma = 3.74;

    std::vector<Point> random_points;
    std::vector<double> random_values;

    std::cout << "\n--- Step 2: Generating random data ---" << std::endl;
    generate_random_data(N, M, sigma, random_points, random_values);
    
    double mean, min_val, max_val, std_dev;
    calculate_statistics(random_values, mean, min_val, max_val, std_dev);

    std::cout << "\nGenerated data statistics:" << std::endl;
    std::cout << "Number of points: " << N << std::endl;
    std::cout << "Theoretical mean: " << M << std::endl;
    std::cout << "Actual mean: " << mean << std::endl;
    std::cout << "Theoretical std: " << sigma << std::endl;
    std::cout << "Actual std: " << std_dev << std::endl;
    std::cout << "Min value: " << min_val << std::endl;
    std::cout << "Max value: " << max_val << std::endl;

    std::cout << "\nUsing ALL " << N << " points for spline analysis" << std::endl;

   
    std::vector<double> p_values = { 0.0, 0.3, 0.7, 0.99 };

    std::cout << "\n--- Step 3: Saving spline data for visualization ---" << std::endl;
    save_spline_comparison(random_points, random_values, p_values);

    analyze_smoothing_parameter(random_points, random_values, p_values);
    return 0;
}
