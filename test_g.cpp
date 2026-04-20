/**
 * @file test_g.cpp
 * @brief Radial distribution function g(r) at four reduced densities.
 *
 * Runs 8 simulations (4 densities × 2 integrators) with rho_star in {0.1, 0.2, 0.3, 0.4}.
 * Both Nosé-Hoover and Brownian dynamics are used to compare structural properties.
 * g(r) data is saved to dati/rdf/NH_0…NH_3.dat and dati/rdf/B_0…B_3.dat.
 */
#include "./sim.hpp"
#include <string>
int main(int argc, char **argv){
  double rho[4]={0.1,0.2,0.3,0.4};

  // NOSE HOOVER
  for(int i=0; i<4; i++){
    mdsim<particleLJ> md;
    md.pars.set_rho_star(rho[i]); // Imposto la densità
    md.prepare_initial_conf(); // Creo la configurazione iniziale (cambia la densità)
    md.pars.set_sim_type(1); // Imposto Nose Hoover
    md.pars.set_rdf(1); // Vogliamo calcolare la g 
    md.pars.set_ssf(0); // Non siamo interessati ad S_k
    std::string filename = "dati/rdf/NH_" + std::to_string(i) + ".dat"; // Nome del file
    md.rad.set_output_name(filename); // Imposto il nuovo nome del file
    md.init_rng(0); // inizializzo il generatore di numeri casuali
    std::cout << "NH:"<<i; 
    md.run(); // lancio la simulazione MD
    std::cout << " OK" << std::endl;

  }

  // DINAMICA BROWNIANA
  for(int i=0; i<4; i++){
    mdsim<particleLJ> md;
    md.pars.set_rho_star(rho[i]); // Imposto la densità
    md.prepare_initial_conf(); // Creo la configurazione iniziale (cambia la densità)
    md.pars.set_sim_type(0); // Imposto dinamica browniana
    md.pars.set_rdf(1); // Vogliamo calcolare la g 
    md.pars.set_ssf(0); // Non siamo interessati ad S_k
    std::string filename = "dati/rdf/B_" + std::to_string(i) + ".dat"; // Nome del file
    md.rad.set_output_name(filename); // Imposto il nuovo nome del file
    md.init_rng(0); // inizializzo il generatore di numeri casuali
    std::cout << "B:"<<i; 
    md.run(); // lancio la simulazione MD
    std::cout << " OK" << std::endl;

  }
  
  return 0;
}
