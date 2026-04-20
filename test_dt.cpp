/**
 * @file test_dt.cpp
 * @brief Time-step convergence test for the Nosé-Hoover integrator.
 *
 * Runs 5 independent simulations with dt in {0.001, 0.00325, 0.0055, 0.00775, 0.01}.
 * Total simulation time is kept fixed; only dt (and thus totsteps) changes.
 * Results are saved to dati/dt/NH_0 … NH_4 for post-processing.
 */
#include "./sim.hpp"
#define N 5

int main(int argc, char **argv){
  double dt[N] = {0.001, 0.00325, 0.0055, 0.00775,  0.01};
  for(int i = 0; i<N; i++){
    mdsim<particleLJ> md;
    md.pars.set_dt(dt[i]);
    md.pars.set_rdf(0);
    md.pars.set_ssf(0); 
    std::string filename = "dati/dt/NH_" + std::to_string(i);
    md.pars.set_output_name(filename);
    md.init_rng(0); 
    md.prepare_initial_conf();
    std::cout << "step: #" << i;
    md.run(); 
    std::cout << " OK" << std::endl;
  }
  return 0;
}
