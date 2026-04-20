/**
 * @file params.hpp
 * @brief Simulation parameters container (`simpars`).
 *
 * Holds all physical and numerical settings for a molecular dynamics run:
 *   - System size: nx * ny * nz particles on a cubic lattice
 *   - Lennard-Jones parameters: sigma, epsilon, cutoff radius rc
 *   - Thermodynamic state: reduced density rho_star, equilibrium temperature T_eq
 *   - Integrator: time step dt, total time t_tot, equilibration steps eqstps
 *   - Thermostat flags: NH (Nosé-Hoover) or B (Brownian dynamics)
 *   - Observable flags: RDF (radial distribution function), SSF (static structure factor)
 *
 * Helper setters (set_dt, set_rho_star, set_sim_type, …) keep derived quantities consistent.
 */
#ifndef _PARAMS_
#define _PARAMS_
#include "./pvector.hpp"
#include <fstream>
#include <string>
#include <cmath>


class simpars{
public:
  using ntype=double;
  int nx, ny, nz; /* nx*ny*nz particelle */
  double T_eq, T_i; // Temperatura di equilibrio, temperatura iniziale
  int Np, Nf; // numero di particelle e numero di gradi di libertà
  long int maxadjstps, eqstps, adjstps;
  long int savemeasure, outstps, totsteps; 
  double rho, rc, rho_star; // density
  int NH, B, RDF, SSF; // Nose-Hover, Brownian dynamics, Radial distribution function, Static structure factor 
  int seed; // -1 means random
  pvector<double, 3> L; // scatola
  double sigma, epsilon, mass, tau; // Parametri Lennard-Jones + tempo tipico
  double dt, t_tot;
  std::string output_name; // Nome del file dove salveremo i dati di output
  simpars(){
    nx = 5;  
    ny = 5;
    nz = 5;
    Np = nx*ny*nz; // Calcoliamo il numero totale di particelle 
    Nf = 3*Np; // Numero di gradi di libertà
    sigma=1;

    // rho_star rappresenta il rapporto tra la densità del sistema e la densità massima.
    // In formule rho_star * rho_max = rho. 
    // Poiché rho_max = (1/sigma)^3 segue che rho = rho_star/(sigma^3)
    rho_star = 0.4; 
    rho = rho_star/pow(sigma, 3);  

    // Variabili che indicano il tipo di simulazione
    NH = 1;
    B = 0;

    // Variabili per indicare se fare la rdf o lo ssf
    RDF = 0;
    SSF = 1;

    // Parametri del potenziale L-J
    epsilon=1;
    rc = 2.5;
    seed=0;
    mass=1;
    
    dt = 0.01; // passo di integrazione
    tau = sigma*sqrt(mass/epsilon); // tempo tipico di un processo diffusivo 
    t_tot = 3*tau;                  // t_tot determinato in funzione di tau  (si sceglie t_tot=3tau)
    totsteps = (long int)rint(t_tot/dt); // Conseguentemente si calcolano gli step totali
    eqstps = (long int)rint(totsteps*(1./3)); // Negativo -> non c'è equilibratura
    
    
    savemeasure=5; 
    outstps=-1;

    T_eq = 1.4; // Temperatura utilizzata nel termostato Nosè-Hoover e nella dinamica browniana
    T_i = 1.4; // Temperatura con cui inizializziamo le velocità iniziali delle particelle

    output_name="dati/results.dat";    // Nome del file che conterrà i dati di output della simulazione
  }

  // Questa funzione permette di cambiare il passo di integrazione, regolando automaticamente il numero di step, per far in modo che t_tot rimanga lo stesso.
  void set_dt(double new_dt){
    dt = new_dt;
    totsteps = (long int)rint(t_tot/new_dt);
    // Se eqstps<0 allora non vogliamo l'equilibratura dentro la nostra simulazione, dunque non la modifichiamo. 
    if (eqstps>0) // se l'equilibratura è minore di 0 non vogliamo modificarla
      eqstps = (long int)rint(totsteps*0.3); 
  }

  // Funzione per cambiare il nome del file di output. Utile per non sovrascrivere files.
  void set_output_name(std::string new_name){
    output_name=new_name;
  }

  // Funzione per cambiare rho_star e contemporaneamente anche la densità rho
  void set_rho_star(double new_rho_star){
    rho_star = new_rho_star;
    rho = rho_star/pow(sigma, 3);  
  }

  // 1 NH, 0 B
  // Con questa funzione cambio il tipo di simulazione in modo tale che non sia possibile che entrambi B e NH siano >=1.
  void set_sim_type(int new_type){
    if (new_type==1){
      NH=1;
      B=0;
    }
    else if (new_type==0){
      NH=0;
      B=1;
    }  
  }

  // Impostiamo se calcolare lo sff
  void set_ssf(int opt){
    SSF = opt;
  }
  
  // Impostiamo se calcolare la rdf
  void set_rdf(int opt){
    RDF = opt;
  }


};
#endif
