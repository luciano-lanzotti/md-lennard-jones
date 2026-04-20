/**
 * @file radial.hpp
 * @brief Radial distribution function g(r) accumulator and evaluator.
 *
 * Usage pattern:
 *   1. Call update_counts(parts, Np, L) at every production step.
 *   2. Call evaluate_g(Np, rho) once at the end of the run.
 *   3. Call save_distribution() to write (r, g(r)) pairs to file.
 *
 * The histogram uses n_bins = 400 bins over [0, r_max=20].
 * Pair distances are computed with minimum-image PBC.
 * Each pair (i,j) contributes 2 to the bin count (accounts for both i→j and j→i).
 */
#ifndef RADIAL_HPP
#define RADIAL_HPP

#include "./pvector.hpp"
#include "./particle.hpp"
#include <cmath>
#include <string>

// Creo una classe per calcolare la distribuzione radiale delle particelle
template<typename particle_type> class radial_distribution{
  long int n_bins; // Numero di bins
  long int n_update; // Numero di volte che viene chiamata la funzione update_counts()
  double r_max; // Distanza massima del centro del sistema
  double dr; // Lunghezza di un bin
  std::vector<long int> counts; // conteggi per ogni bin
  std::vector<double> g; // Radial distribution function 
  std::string output_name; // Nome del file su cui salveremo la g(r)
public:
  radial_distribution(){
    output_name="dati/radial.dat";
    n_update = 0; // All'inizio va sempre posto n_updato pari a 0
    r_max = 20.; // Si osserva che per r maggiori di 20 la g(r) è sempre nulla 
    n_bins = 400; 
    dr = r_max/n_bins;
    counts.resize(n_bins); // Array che conterrà il numero di particelle con distanza in (r,r+dr)
    g.resize(n_bins); // Array che conterrà la g(r)
    for(int i=0; i<n_bins; i++)
      counts[i]=0; //Serve inizializzare a 0 i conteggi 
  }

  // Passo per riferimento per non creare ogni volta una copia del vettore particelle.
  // update_counts() è pensata per essere chiamata ad ogni step temporale (passata l'equlibratura).
  void update_counts(std::vector<particle_type>& parts, long int Np, pvector<double,3> L){
    double rn; // distanza tra particelle
    long int bin; //bin da incrementare
    pvector<double,3> Dr; 
    // Ciclo for per j>i
    for (int i=0; i < Np; i++){
      for (int j=i+1; j < Np; j++){
        Dr = parts[i].r - parts[j].r; // Calcolo della distanza 
        Dr -= L.mulcw(rint(Dr.divcw(L)));// Applicazione delle PBC
        rn = Dr.norm();
        // Solo rn<r_max aggiorniamo l'elemento count[bin]
        // Se rn>r_max la distanza è fuori l'intervallo considerato
        if (rn<r_max){
          bin = (long int)(rint((rn/r_max)*n_bins));
          counts[bin] += 2;
        }
      }
    }
    n_update += 1; // Poiché questa funzione viene chiamata nel run() solo dopo l'equilibratura, 
    // a priori non sappiamo quante volte aggiorneremo l'istogramma.
    // Quindi non sappiamo per quale valore dividere l'array counts per ottenere il numero medio di particelle in un dato bin.
  }

  // evaluate_g() è pensata per essere chiamata alla fine della simulazione
  void evaluate_g(long int Np, double rho){
    double average_count; // Rappresenta il numero medio di particelle in un dato bin
    double V_shell; // Volume del guscio sferico infinitesimo
    for(int i=0; i<n_bins; i++){ // Ciclo for su tutti i bin
      average_count = (double)counts[i]/n_update;
      V_shell = 4*M_PI*dr*pow((i+0.5)*dr,2); // 4*pi*(r^2)*dr. Si sceglie r come il valore al centro del bin, cioè (i+0.5)*dr
      g[i] = average_count/(V_shell*rho*Np); // Formula della g(r)
    }
  }
  
  // save_distribution() è pensata per essere chiamata alla fine delle simulazione
  void save_distribution(){
    std::fstream f;
    f.open(output_name,std::ios::out|std::ios::trunc); 
    for(int i=0; i<n_bins; i++){
      f << dr*(i+0.5) << "\t "; // Salviamo r come il valore centrale del bin
      f << g[i] << "\n ";
    }
    f.sync();      
    f.close();
  }

  // Funzione per cambiare il nome del file di output
  // Questo metodo torna utile quando si lanciano più simulazioni e non si vogliono sovrascrivere i risultati
  void set_output_name(std::string new_name){
    output_name = new_name;
  }

};
#endif