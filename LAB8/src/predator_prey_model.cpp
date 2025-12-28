#include "../include/predator_prey_model.h"
#include <fstream>
#include <iostream>
#include <iomanip>

population_dynamics::population_dynamics(const model_parameters& p) : params(p) {}

void population_dynamics::compute_derivatives(double prey_count, double predator_count,
                                             double& dprey_dt, double& dpredator_dt) const {
    dprey_dt = (params.prey_birth_rate - params.hunting_success_rate * predator_count) * prey_count;
    dpredator_dt = (-params.predator_death_rate + params.predator_reproduction_rate * prey_count) * predator_count;
}

void population_dynamics::get_equilibrium_point(double& equilibrium_prey, 
                                               double& equilibrium_predator) const {
    equilibrium_prey = params.predator_death_rate / params.predator_reproduction_rate;
    equilibrium_predator = params.prey_birth_rate / params.hunting_success_rate;
}

void population_dynamics::save_to_file(const std::string& filename, double* prey_data,
                                      double* predator_data, double* time_data, int size) const {
    std::ofstream output_file(filename);
    if (!output_file.is_open()) {
        std::cerr << "Ошибка открытия файла: " << filename << std::endl;
        return;
    }
    
    output_file << std::fixed << std::setprecision(6);
    output_file << "t,x,y\n";
    for (int i = 0; i < size; ++i) {
        output_file << time_data[i] << "," << prey_data[i] << "," << predator_data[i] << "\n";
    }
    
    output_file.close();
}