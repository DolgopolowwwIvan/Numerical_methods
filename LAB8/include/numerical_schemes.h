#ifndef NUMERICAL_SCHEMES_H
#define NUMERICAL_SCHEMES_H

#include "predator_prey_model.h"
#include <vector>

class time_integrator {
protected:
    double time_step;
    int total_steps;
    
public:
    time_integrator(double dt, int steps);
    virtual ~time_integrator() = default;
    
    virtual void integrate(const population_dynamics& model,
                          double initial_prey, double initial_predator,
                          std::vector<double>& prey_history,
                          std::vector<double>& predator_history,
                          std::vector<double>& time_history) = 0;
};

class runge_kutta_4th_order : public time_integrator {
public:
    runge_kutta_4th_order(double dt, int steps);
    
    void integrate(const population_dynamics& model,
                  double initial_prey, double initial_predator,
                  std::vector<double>& prey_history,
                  std::vector<double>& predator_history,
                  std::vector<double>& time_history) override;
};

class adams_bashforth_3rd_order : public time_integrator {
public:
    adams_bashforth_3rd_order(double dt, int steps);
    
    void integrate(const population_dynamics& model,
                  double initial_prey, double initial_predator,
                  std::vector<double>& prey_history,
                  std::vector<double>& predator_history,
                  std::vector<double>& time_history) override;
};

#endif 