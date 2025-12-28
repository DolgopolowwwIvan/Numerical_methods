#include "../include/numerical_schemes.h"
#include <vector>
#include <fstream>
#include <iostream>

time_integrator::time_integrator(double dt, int steps) 
    : time_step(dt), total_steps(steps) {}

runge_kutta_4th_order::runge_kutta_4th_order(double dt, int steps) 
    : time_integrator(dt, steps) {}

void runge_kutta_4th_order::integrate(const population_dynamics& model,
                                     double initial_prey, double initial_predator,
                                     std::vector<double>& prey_history,
                                     std::vector<double>& predator_history,
                                     std::vector<double>& time_history) {
    
    prey_history.clear();
    predator_history.clear();
    time_history.clear();
    
    prey_history.reserve(total_steps + 1);
    predator_history.reserve(total_steps + 1);
    time_history.reserve(total_steps + 1);
    
    double current_prey = initial_prey;
    double current_predator = initial_predator;
    double current_time = 0.0;
    
    prey_history.push_back(current_prey);
    predator_history.push_back(current_predator);
    time_history.push_back(current_time);
    
    for (int step = 0; step < total_steps; ++step) {
        double k1_prey, k1_predator;
        double k2_prey, k2_predator;
        double k3_prey, k3_predator;
        double k4_prey, k4_predator;
        
        model.compute_derivatives(current_prey, current_predator, k1_prey, k1_predator);
        model.compute_derivatives(current_prey + time_step * k1_prey / 2.0,
                                 current_predator + time_step * k1_predator / 2.0,
                                 k2_prey, k2_predator);
        model.compute_derivatives(current_prey + time_step * k2_prey / 2.0,
                                 current_predator + time_step * k2_predator / 2.0,
                                 k3_prey, k3_predator);
        model.compute_derivatives(current_prey + time_step * k3_prey,
                                 current_predator + time_step * k3_predator,
                                 k4_prey, k4_predator);
        
        current_prey += time_step * (k1_prey + 2.0 * k2_prey + 2.0 * k3_prey + k4_prey) / 6.0;
        current_predator += time_step * (k1_predator + 2.0 * k2_predator + 2.0 * k3_predator + k4_predator) / 6.0;
        current_time += time_step;
        
        prey_history.push_back(current_prey);
        predator_history.push_back(current_predator);
        time_history.push_back(current_time);
    }
}

adams_bashforth_3rd_order::adams_bashforth_3rd_order(double dt, int steps) 
    : time_integrator(dt, steps) {}

void adams_bashforth_3rd_order::integrate(const population_dynamics& model,
                                         double initial_prey, double initial_predator,
                                         std::vector<double>& prey_history,
                                         std::vector<double>& predator_history,
                                         std::vector<double>& time_history) {
    
    prey_history.clear();
    predator_history.clear();
    time_history.clear();
    
    const int history_size = 3; 
 
    std::vector<double> x_history(history_size);
    std::vector<double> y_history(history_size);
    std::vector<double> fx_history(history_size);
    std::vector<double> fy_history(history_size);
    
    double x = initial_prey;
    double y = initial_predator;
    double t = 0.0;
    
 
    prey_history.push_back(x);
    predator_history.push_back(y);
    time_history.push_back(t);
    
    double fx, fy;
    model.compute_derivatives(x, y, fx, fy);
    x_history[0] = x;
    y_history[0] = y;
    fx_history[0] = fx;
    fy_history[0] = fy;

    for (int i = 1; i < history_size; ++i) {
        double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;
        
    
        model.compute_derivatives(x, y, k1x, k1y);
        
     
        model.compute_derivatives(x + time_step * k1x / 2.0,
                                 y + time_step * k1y / 2.0,
                                 k2x, k2y);
        
     
        model.compute_derivatives(x + time_step * k2x / 2.0,
                                 y + time_step * k2y / 2.0,
                                 k3x, k3y);
        
       
        model.compute_derivatives(x + time_step * k3x,
                                 y + time_step * k3y,
                                 k4x, k4y);
        
       
        x += time_step * (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0;
        y += time_step * (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0;
        t += time_step;
        
        prey_history.push_back(x);
        predator_history.push_back(y);
        time_history.push_back(t);

        if (i < history_size) {
            x_history[i] = x;
            y_history[i] = y;
            model.compute_derivatives(x, y, fx_history[i], fy_history[i]);
        }
    }
    
    // Основной цикл Adams-Bashforth 3-го порядка
    for (int step = history_size - 1; step < total_steps; ++step) {

        double x_new = x + time_step / 12.0 * 
            (23.0 * fx_history[2] - 16.0 * fx_history[1] + 5.0 * fx_history[0]);
        
        double y_new = y + time_step / 12.0 * 
            (23.0 * fy_history[2] - 16.0 * fy_history[1] + 5.0 * fy_history[0]);
        
        t += time_step;
        

        prey_history.push_back(x_new);
        predator_history.push_back(y_new);
        time_history.push_back(t);

        for (int i = 0; i < history_size - 1; ++i) {
            x_history[i] = x_history[i + 1];
            y_history[i] = y_history[i + 1];
            fx_history[i] = fx_history[i + 1];
            fy_history[i] = fy_history[i + 1];
        }

        x = x_new;
        y = y_new;
        x_history[2] = x;
        y_history[2] = y;
        
        model.compute_derivatives(x, y, fx_history[2], fy_history[2]);
    }
}