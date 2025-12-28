#ifndef PREDATOR_PREY_MODEL_H
#define PREDATOR_PREY_MODEL_H

#include <string>

struct model_parameters {
    double prey_birth_rate;       // коэффициент рождаемости жертв
    double hunting_success_rate;  // коэффициент успешной охоты
    double predator_death_rate;   // коэффициент убыли хищников
    double predator_reproduction_rate;   // коэффициент воспроизводства хищников
};

class population_dynamics {
private:
    model_parameters params;
    
public:
    population_dynamics(const model_parameters& p);
    
    void compute_derivatives(double prey_count, double predator_count, 
                            double& dprey_dt, double& dpredator_dt) const;
    
    void get_equilibrium_point(double& equilibrium_prey, double& equilibrium_predator) const;
    
    void save_to_file(const std::string& filename, double* prey_data, 
                     double* predator_data, double* time_data, int size) const;
};

#endif 