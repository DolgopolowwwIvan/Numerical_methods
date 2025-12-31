#include "../include/predator_prey_model.h"
#include "../include/numerical_schemes.h"
#include <iostream>
#include <vector>
#include <memory>

int main() {
    
    model_parameters lions_antelopes = {
        0.007,      // α - в 10 раз меньше
        0.00035,    // β - в 10 раз меньше
        0.00174,    // γ - в 10 раз меньше  
        0.000058    // δ - в 10 раз меньше
    };

    model_parameters dragonfly_mosquito = {
        0.8,      // α
        0.02,     // β
        0.6,      // γ
        0.01      // δ
    };
    
    population_dynamics savanna_system(lions_antelopes);
    population_dynamics aquatic_system(dragonfly_mosquito);
    
    double antelopes_eq, lions_eq;
    savanna_system.get_equilibrium_point(antelopes_eq, lions_eq);
    
    double mosquitoes_eq, dragonflies_eq;
    aquatic_system.get_equilibrium_point(mosquitoes_eq, dragonflies_eq);
    // Около равновесия
    double antelopes_case1 = 32.0;
    double lions_case1 = 19.0;
    
    //  Больше хищников 
    double antelopes_case2 = 11.0;
    double lions_case2 = 29.0;
    
    // Больше жертв 
    double antelopes_case3 = 31.0;
    double lions_case3 = 2.0;


    //Много жертв
    double mosquitoes_case1 = 55.0;    
    double dragonflies_case1 = 15.0;     

    //Много хищников
    double mosquitoes_case2 = 6.0;     
    double dragonflies_case2 = 15.0;    
    
    //Близко к равновесию
    double mosquitoes_case3 = 66.0;     
    double dragonflies_case3 = 44.0;    
    

    double dt_savanna = 1.0;
    int steps_savanna = 3650;

    double dt_aquatic = 0.1;
    int steps_aquatic = 1600;  // 160 дней × 10 шагов в день

    auto rk4_savanna = std::make_unique<runge_kutta_4th_order>(dt_savanna, steps_savanna);
    auto ab3_savanna = std::make_unique<adams_bashforth_3rd_order>(dt_savanna, steps_savanna);
    
    auto rk4_aquatic = std::make_unique<runge_kutta_4th_order>(dt_aquatic, steps_aquatic);
    auto ab3_aquatic = std::make_unique<adams_bashforth_3rd_order>(dt_aquatic, steps_aquatic);
    
    std::vector<double> prey_rk4, predator_rk4, time_rk4;
    std::vector<double> prey_ab3, predator_ab3, time_ab3;
    

    rk4_savanna->integrate(savanna_system, antelopes_case1, lions_case1,
                          prey_rk4, predator_rk4, time_rk4);
    ab3_savanna->integrate(savanna_system, antelopes_case1, lions_case1,
                          prey_ab3, predator_ab3, time_ab3);
    
    savanna_system.save_to_file("lions_antelopes_case1_rk4.csv", 
                               prey_rk4.data(), predator_rk4.data(), 
                               time_rk4.data(), prey_rk4.size());
    savanna_system.save_to_file("lions_antelopes_case1_ab3.csv",
                               prey_ab3.data(), predator_ab3.data(),
                               time_ab3.data(), prey_ab3.size());
    

    prey_rk4.clear(); predator_rk4.clear(); time_rk4.clear();
    prey_ab3.clear(); predator_ab3.clear(); time_ab3.clear();
    
    rk4_savanna->integrate(savanna_system, antelopes_case2, lions_case2,
                          prey_rk4, predator_rk4, time_rk4);
    ab3_savanna->integrate(savanna_system, antelopes_case2, lions_case2,
                          prey_ab3, predator_ab3, time_ab3);
    
    savanna_system.save_to_file("lions_antelopes_case2_rk4.csv",
                               prey_rk4.data(), predator_rk4.data(),
                               time_rk4.data(), prey_rk4.size());
    savanna_system.save_to_file("lions_antelopes_case2_ab3.csv",
                               prey_ab3.data(), predator_ab3.data(),
                               time_ab3.data(), prey_ab3.size());
    

    prey_rk4.clear(); predator_rk4.clear(); time_rk4.clear();
    prey_ab3.clear(); predator_ab3.clear(); time_ab3.clear();
    
    rk4_savanna->integrate(savanna_system, antelopes_case3, lions_case3,
                          prey_rk4, predator_rk4, time_rk4);
    ab3_savanna->integrate(savanna_system, antelopes_case3, lions_case3,
                          prey_ab3, predator_ab3, time_ab3);
    
    savanna_system.save_to_file("lions_antelopes_case3_rk4.csv",
                               prey_rk4.data(), predator_rk4.data(),
                               time_rk4.data(), prey_rk4.size());
    savanna_system.save_to_file("lions_antelopes_case3_ab3.csv",
                               prey_ab3.data(), predator_ab3.data(),
                               time_ab3.data(), prey_ab3.size());

    prey_rk4.clear(); predator_rk4.clear(); time_rk4.clear();
    prey_ab3.clear(); predator_ab3.clear(); time_ab3.clear();
    
    rk4_aquatic->integrate(aquatic_system, mosquitoes_case1, dragonflies_case1,
                          prey_rk4, predator_rk4, time_rk4);
    ab3_aquatic->integrate(aquatic_system, mosquitoes_case1, dragonflies_case1,
                          prey_ab3, predator_ab3, time_ab3);
    
    aquatic_system.save_to_file("dragonfly_mosquito_case1_rk4.csv",
                               prey_rk4.data(), predator_rk4.data(),
                               time_rk4.data(), prey_rk4.size());
    aquatic_system.save_to_file("dragonfly_mosquito_case1_ab3.csv",
                               prey_ab3.data(), predator_ab3.data(),
                               time_ab3.data(), prey_ab3.size());
    

    prey_rk4.clear(); predator_rk4.clear(); time_rk4.clear();
    prey_ab3.clear(); predator_ab3.clear(); time_ab3.clear();
    
    rk4_aquatic->integrate(aquatic_system, mosquitoes_case2, dragonflies_case2,
                          prey_rk4, predator_rk4, time_rk4);
    ab3_aquatic->integrate(aquatic_system, mosquitoes_case2, dragonflies_case2,
                          prey_ab3, predator_ab3, time_ab3);
    
    aquatic_system.save_to_file("dragonfly_mosquito_case2_rk4.csv",
                               prey_rk4.data(), predator_rk4.data(),
                               time_rk4.data(), prey_rk4.size());
    aquatic_system.save_to_file("dragonfly_mosquito_case2_ab3.csv",
                               prey_ab3.data(), predator_ab3.data(),
                               time_ab3.data(), prey_ab3.size());
    

    prey_rk4.clear(); predator_rk4.clear(); time_rk4.clear();
    prey_ab3.clear(); predator_ab3.clear(); time_ab3.clear();
    
    rk4_aquatic->integrate(aquatic_system, mosquitoes_case3, dragonflies_case3,
                          prey_rk4, predator_rk4, time_rk4);
    ab3_aquatic->integrate(aquatic_system, mosquitoes_case3, dragonflies_case3,
                          prey_ab3, predator_ab3, time_ab3);
    
    aquatic_system.save_to_file("dragonfly_mosquito_case3_rk4.csv",
                               prey_rk4.data(), predator_rk4.data(),
                               time_rk4.data(), prey_rk4.size());
    aquatic_system.save_to_file("dragonfly_mosquito_case3_ab3.csv",
                               prey_ab3.data(), predator_ab3.data(),
                               time_ab3.data(), prey_ab3.size());
    
    return 0;
}