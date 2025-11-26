#include <cmath>
#include "integral_calculators.h"

double TrapezoidalRule::compute(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i) sum += f(a + i * h);
    return sum * h;
}

double SimpsonRule::compute(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += (i % 2 == 0 ? 2 : 4) * f(x);
    }
    return sum * h / 3.0;
}

double RichardsonExtrapolation::compute(double Ih, double Ih2, int k) {
    return Ih2 + (Ih2 - Ih) / (std::pow(2, k) - 1);
}