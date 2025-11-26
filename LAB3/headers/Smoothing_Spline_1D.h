#pragma once
#ifndef SMOOTHING_SPLINE_1D_H
#define SMOOTHING_SPLINE_1D_H

#include "Spline.h"

namespace Com_Methods {
    class Smoothing_Spline_1D : public Spline {
    private:
        double SMOOTH;  // параметр сглаживания p
        std::vector<Point> Points;
        std::vector<double> alpha; // коэффициенты разложения
        
        void Transition_To_Master_Element(int Seg_Num, const double &X, double &Ksi) const;
        double Basis_Function(int Number, const double &Ksi) const;
        double Der_Basis_Function(int Number, const double &Ksi) const;
        
    public:
        Smoothing_Spline_1D(const double &SMOOTH = 0.5);
        void Update_Spline(const std::vector<Point> &Points, 
                          const std::vector<double> &F_Value) override;
        void Get_Value(const Point &P, double *Res) const override;
        
        // Дополнительные методы для удобства
        void Set_Smoothing_Parameter(double p) { SMOOTH = p; }
        double Get_Smoothing_Parameter() const { return SMOOTH; }
    };
}

#endif