/**
 * @file test_U.cpp
 * @brief Energy conservation test for both Nosé-Hoover and Brownian dynamics.
 *
 * Sweeps 5 time steps in {0.005, 0.01125, 0.0175, 0.02375, 0.03} for each integrator.
 * Tracks the conserved energy (K + U_shifted + E_thermostat) to evaluate integrator accuracy.
 * Results saved to dati/u/NH_0…NH_4 and dati/u/B_0…B_4.
 */
#include "./sim.hpp"
#define N 5

int main(int argc, char **argv){
  double dt[N] = {0.005  , 0.01125, 0.0175 , 0.02375, 0.03   };

  // NH
  for(int i = 0; i<N; i++){
    mdsim<particleLJ> md;
    md.pars.set_sim_type(1);
    md.pars.set_dt(dt[i]);
    md.pars.set_rdf(0);
    md.pars.set_ssf(0);
    std::string filename = "dati/u/NH_" + std::to_string(i);
    md.pars.set_output_name(filename);
    md.init_rng(0); 
    md.prepare_initial_conf();
    std::cout << "step (NH): #" << i;
    md.run(); 
    std::cout << " OK" << std::endl;
  }

  // Brownian 
  for(int i = 0; i<N; i++){
    mdsim<particleLJ> md;
    md.pars.set_sim_type(0);
    md.pars.set_dt(dt[i]);
    md.pars.set_rdf(0);
    md.pars.set_ssf(0);
    std::string filename = "dati/u/B_" + std::to_string(i);
    md.pars.set_output_name(filename);
    md.init_rng(0); 
    md.prepare_initial_conf();
    std::cout << "step (B): #" << i;
    md.run(); 
    std::cout << " OK" << std::endl;
  }
  return 0;
}
