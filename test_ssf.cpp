/**
 * @file test_ssf.cpp
 * @brief Static structure factor S(k) at four higher reduced densities.
 *
 * Runs 8 simulations (4 densities × 2 integrators) with rho_star in {0.5, 0.6, 0.7, 0.8}.
 * S(k) is computed by randomly sampling wave vectors k = (2pi/L)*n and averaging |rho(k)|^2.
 * Results saved to dati/ssf/NH_0…NH_3.dat and dati/ssf/B_0…B_3.dat.
 */
#include "./sim.hpp"
#include <string>
int main(int argc, char **argv){
  double rho[4]={0.5,0.6,0.7,0.8};

  // NOSE HOOVER
  for(int i=0; i<4; i++){
    mdsim<particleLJ> md;
    md.pars.set_rho_star(rho[i]);
    md.prepare_initial_conf(); // Creo la configurazione iniziale (cambia la densità)
    md.pars.set_sim_type(1); // Imposto Nose Hoover
    md.pars.set_rdf(0); // NON calcoliamo la radial distribution function
    md.pars.set_ssf(1); // Calcoliamo lo static structure factor
    std::string filename = "dati/ssf/NH_" + std::to_string(i) + ".dat"; // Nome del file
    md.ssf.set_output_name(filename); // impostiamo il nome del file (ci saranno due file per ogni densita', uno per NH e l'altro per la dinamica browniana)
    md.init_rng(0); // inizializzo il generatore di numeri casuali
    std::cout << "NH:"<<i; 
    md.run(); // lancio la simulazione MD
    std::cout << " OK" << std::endl;
  }

  // DINAMICA BROWNIANA
  for(int i=0; i<4; i++){
    mdsim<particleLJ> md;
    md.pars.set_rho_star(rho[i]);
    md.prepare_initial_conf(); // Creo la configurazione iniziale (cambia la densità)
    md.pars.set_sim_type(0); // Imposto dinamica browniana
    md.pars.set_rdf(0); // NON calcoliamo la radial distribution function
    md.pars.set_ssf(1); // calcoliamo lo static structure factor
    std::string filename = "dati/ssf/B_" + std::to_string(i) + ".dat"; // Nome del file
    md.ssf.set_output_name(filename); // impostiamo il nome del file (ci saranno due file per ogni densita', uno per NH e l'altro per la dinamica browniana)
    md.init_rng(0); // inizializzo il generatore di numeri casuali
    std::cout << "B:"<<i; 
    md.run(); // lancio la simulazione MD
    std::cout << " OK" << std::endl;
  }
  
  return 0;
}
