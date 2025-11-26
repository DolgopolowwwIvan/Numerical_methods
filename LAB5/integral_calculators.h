#ifndef INTEGRAL_CALCULATORS_H
#define INTEGRAL_CALCULATORS_H

class TrapezoidalRule {
public:
    double compute(double (*f)(double), double a, double b, int n);
    int getOrder() const { return 2; }
};

class SimpsonRule {
public:
    double compute(double (*f)(double), double a, double b, int n);
    int getOrder() const { return 4; }
};

class RichardsonExtrapolation {
public:
    static double compute(double Ih, double Ih2, int k);
};

#endif