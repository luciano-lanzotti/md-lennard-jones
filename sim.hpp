/**
 * @file sim.hpp
 * @brief Core MD simulation engine: base `sim` and derived `mdsim` classes.
 *
 * `sim<particle_type>` (base):
 *   - Holds the particle array and simpars parameters.
 *   - prepare_initial_conf(): places particles on a cubic lattice scaled to rho.
 *   - pbc(i): applies periodic boundary conditions to particle i.
 *   - init_rng(seed): seeds the global RNG.
 *
 * `mdsim<particle_type>` (derived from sim):
 *   - prepare_initial_conf(): extends the base with Maxwell-Boltzmann velocity
 *     initialisation (zero net momentum) and SSF setup.
 *   - run(): main loop supporting three integrators selected by pars.NH / pars.B:
 *       - NH=1 : velocity-Verlet + Nosé-Hoover thermostat (Trotter factored)
 *       - B=1  : stochastic Brownian dynamics (Langevin)
 *       - else : plain velocity-Verlet (microcanonical NVE)
 *   - After equilibration, accumulates g(r) and/or S(k) if enabled.
 *   - Output is written to the file specified by pars.output_name.
 */
#ifndef _SIMCLASS_
#define _SIMCLASS_

#include <vector>
#include <string>
#include <fstream>
#include "./params.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include "./nh.hpp"
#include "./brownian.hpp"
#include "./radial.hpp"
#include "./static.hpp"
#include <iomanip> // for setprecision()
#include <math.h>
#define M_PI 3.14159265358979323846
using ntype=double;

template<typename particle_type>
class sim{
  
protected: 
  std::vector<particle_type> parts; // Vettore di particelle

  void pbc(int i){
    auto Dr = parts[i].r;
    Dr = pars.L.mulcw(rint(Dr.divcw(pars.L))); // L*rint(Dr/L)
    parts[i].r = parts[i].r - Dr;
  }
  

public:
  using simp = simpars;
  simp pars;  // Parametri generali della simulazione

  void prepare_initial_conf(void){
    int ix, iy, iz;
    int cc=0;
    parts.resize(pars.Np);
    ntype vcell = pow(pars.sigma,3.0);
    ntype rhomax = 1.0/vcell;
    ntype sf;
    sf = cbrt(rhomax/pars.rho);
    pars.L ={ntype(pars.nx), ntype(pars.ny), ntype(pars.nz)};
    ntype clen = sf*pars.sigma;
    pars.L *= clen;
    for (ix = 0; ix < pars.nx; ix++)
      for (iy = 0; iy < pars.ny; iy++)
        for (iz = 0; iz < pars.nz; iz++)
          {
            parts[cc].r = {ix*clen, iy*clen, iz*clen};
            parts[cc].set_sigma(pars.sigma);
            parts[cc].set_epsilon(pars.epsilon);
            parts[cc].set_rcut(pars.rc);
            cc++;
            
          }
  }
  
  void init_rng(int n){
    if (n < 0)
      rng.rseed();
    else
      rng.seed(n);
  }

  void run(void){}; 
};


template<typename particle_type> class mdsim: public sim<particle_type>{
  using ntype=double;
  using bc = sim<particle_type>;
  using bc::parts, bc::pbc;
  ntype Us=0; // Energia potenziale riscalata
  ntype U=0; // Energia potenziale
  nh<particle_type> thermostat; // Termostato secondo NH
  brownian<particle_type> b; // Classe che implementa l'evolutore exp^{iL_G}
  
  // La funzione per inizializzare le misure scrive una prima riga dove riassume i parametri importanti della simulazione.
  // Con l'asterisco python numpy ignorerà la prima riga.
  void init_measures(void){
    std::fstream f;
    f.open(pars.output_name, std::ios::out|std::ios::trunc);
    f << "#dt:" << pars.dt;
    f << " totstep:" << pars.totsteps;
    f << " tau:" << pars.tau;
    f << " t_tot:" << pars.t_tot;
    f << " rho:" << pars.rho;
    f << " NH:" << pars.NH;
    f << " B:" << pars.B;
    f << std::endl;
    f.close();
  }

  // Calcolo dell'energia cinetica
  ntype calcK(void){
    ntype K=0.0;
    for (int i=0; i < pars.Np; i++){
      K += 0.5*parts[i].m*parts[i].v*parts[i].v;
    }
    return K;
  }

  // Calcolo della temperatura a partire dall'energia cinetica. 
  ntype calcT(void){
    return calcK()*2/pars.Nf; // pars.Nf sono i gradi di libertà del sistema, posti uguali a 3N.
  }
 
  // Funzione per salvare le misure nel formato:
  // t  i*dt  K+Us  E_NH  T  
  void save_measures(long int t){
    std::fstream f;
    f.open(pars.output_name, std::ios::out|std::ios::app);
    ntype K=calcK(); // Energia cinetica
    double E_eta = thermostat.energy(); //P_eta^2/(2*Q) + 3*N*k_BT*eta
    f << t << "\t"; // Passo della simulazione
    f << std::fixed << std::setprecision(3)  << t*pars.dt << "\t"; // tempo della simulazione 
    f << std::setprecision(15) << K+Us << "\t"; // Energia interna
    f << std::setprecision(15) << K+Us+E_eta<< "\t"; // Energia conservata con NH
    f << std::setprecision(15) << calcT() << std::endl;
    f.sync();      
    f.close();
  }

  
  void calc_forces(ntype& Us, ntype& U){
    // Azzero il potenziale per ricalcolarlo da capo
    U  = 0;
    Us = 0; 
    // Azzero anche le forze su ogni particella 
    for (int i=0; i<pars.Np; i++){
      parts[i].f={0,0,0};
    }
    pvector<ntype,3> fij_temp;
    // Due cicli for annidati tali che j>i
    for(int i=0; i<pars.Np; i++){
      for(int j=i+1; j<pars.Np; j++){
        fij_temp = parts[i].fij(parts[j], pars.L); // Forza tra i e j
        // La forza tra i e j sarà positiva se agisce su i, ma negativa se agisce su j.
        parts[i].f += fij_temp; 
        parts[j].f -= fij_temp;
        // Il potenziale è la somma dei contributi dalle interazioni tra i e j
        U  += parts[i].vij(parts[j],pars.L);
        Us += parts[i].vijs(parts[j],pars.L);  
      }
    }


  }

  // Generatore gaussiano con algoritmo box muller 
  ntype gauss(void){
    double x1, x2;
    do
      {
        x1 = rng.ranf();
        x2 = rng.ranf();
      }
    while (x1==0 || x2==0);
    return cos(2*M_PI*x2)*sqrt(-2.0*log(x1));
  }

  // Attraverso questo metodo implementiamo un passo del velocity verlet su tutte le particelle.
  // Avere un metodo così tornerà utile per fare l'algoritmo NH
  void velocity_verlet(void){
    for(int i=0; i<pars.Np; i++){
      parts[i].expiLp(pars.dt*0.5);
      parts[i].expiLq(pars.dt);
      pbc(i);
    }
    calc_forces(Us,U);
    for(int i=0; i<pars.Np; i++){
      parts[i].expiLp(pars.dt*0.5);
    }
  }

public:
using bc::pars; // I parametri della simulazione sono in public
radial_distribution<particle_type> rad; // Classe che calcola la distribuzione radiale
ssf<particle_type> ssf; // Clase che calcola lo static structure factor

  void prepare_initial_conf(void){
    bc::prepare_initial_conf();
    pvector<double,3> total_momentum = {0.0, 0.0, 0.0}; // Creaiamo un vettore momento totale inizializzato a 0
    ntype M=0; // Questa è la massa totale del sistema, la inizializziamo a 0.
    // Ciclo su tutte le particelle
    for (int i=0; i<pars.Np; i++){
      M+=parts[i].m; // Incrementiamo la massa del sistema
      // Ciclo sulle dimensioni spaziali 
      for(int x=0; x<3; x++){ 
        parts[i].v(x) = gauss() * sqrt(pars.T_i/parts[i].m); // Assegnami la velocità alla i-esima particella su una delle 3 direzioni spaziali
        total_momentum(x) += parts[i].v(x)*parts[i].m; // Sommiamo l'impulso appena generato al momento totale
      }
      parts[i].set_vcut(); // Calcoliamo il potenziale di cutoff sulla particella i-esima. 
      // Infatti vcut è un attributo della classe particle e deve essere specificato per tutte le particelle
    }
    // Così otteniamo la velocità del centro di massa. Usiamo la stessa variabile per semplicità.
    total_momentum/=M; 
    // Riscaliamo le velocità di tutte le particelle per la velocità del centro di massa
    for(int i=0; i<pars.Np; i++){
      parts[i].v -= total_momentum;
    }     

    // Prepariamo lo static structure factor
    // Infatti lo ssf per definire alcuni suoi elementi ha bisogno dei valori dei lati della scatola. 
    // Valori disponibili solo dopo aver fatto i passaggi di sopra.
    ssf.init_ssf(pars.L);

    
  }


  void run(void){
    calc_forces(Us, U); // Calcoliamo forze e potenziale per la configurazione iniziale
    init_measures();
    save_measures(0); // Salviamo le misurazioni al tempo zero
    // La variabile t parte da 1 poiché idealmente il tempo 0 è quello iniziale
    for (long int t=1; t < pars.totsteps+1; t++){ 
      // Se abbiamo indicato pars.NH >= allora si farà un passo NH
      if(pars.NH){
        // Gli evolutori della classe nh evolvono tutte le particelle in una sola chiamata della funzione
        thermostat.expiLpeta(pars.dt*0.5, parts);
        thermostat.expiLFNH(pars.dt*0.5, parts);
        thermostat.expiLeta(pars.dt*0.5);        
        velocity_verlet();    
        thermostat.expiLeta(pars.dt*0.5);
        thermostat.expiLFNH(pars.dt*0.5, parts);
        thermostat.expiLpeta(pars.dt*0.5, parts);
      }

      // Altrimenti se abbiamo indicato pars.B >= 1 si farà la dinamica browniana
      else if(pars.B){
        for(int i=0; i<pars.Np; i++){
          parts[i].expiLp(pars.dt*0.5); // Evolutore degli impulsi 
          parts[i].expiLq(pars.dt*0.5); // Evolutore delle posizioni
          pbc(i); // Dopo aver cambiato le posizioni è necessario applicare le pbc
        }
        calc_forces(Us, U); // Dopo aver aggiornato le posizioni le forze sulle singolo particelle cambiano
        for(int i=0; i<pars.Np; i++){
          b.expiLG(pars.dt, parts[i]); // Evolutore proprio della dinamica browniana
          parts[i].expiLq(pars.dt*0.5); // Evolutore delle posizioni
          pbc(i); // Come prima si applicano le pbc
        }
        calc_forces(Us, U); // Ricalcoliamo le forze
        for(int i=0; i<pars.Np; i++)
          parts[i].expiLp(pars.dt*0.5); // Aggiornamento degli impulsi
      }

      // Se entrambi pars.B e pars.NH sono 0, allora evolviamo il sistema con verlet
      else velocity_verlet();



      // L'aggiornamento delle variabili per fare la RDF o lo SSF deve avvenire solo dopo l'equilibratura.
      if(t>pars.eqstps){
        if(pars.RDF) rad.update_counts(parts, pars.Np, pars.L);
        if(pars.SSF) ssf.update_ssf(parts, pars.Np, pars.L);
      }
      
      

      // Stampa su terminale
      if (pars.outstps > 0 && t % pars.outstps== 0){
        ntype K=calcK();
        std::cout << "step #" << t << " K:";
        std::cout << K;
        
        // Questo if permette di visualizzare se il centro di massa è fermo, come dovrebbe essere, o sta traslando erroneamente 
        # if 1
        pvector<ntype,3> momento_totale={0,0,0};
        for(int i=0; i<pars.Np;i++){
          momento_totale += parts[i].v * parts[i].m;
        }
        std::cout << std::setprecision(8) << " P:" << momento_totale;
        #endif

        // Questo if permette di visualizzare i valori delle variabili del termostato NH
        #if 1
        std::cout << " eta:" << std::setprecision(8) << thermostat.eta;
        std::cout << " P_eta:" << std::setprecision(8) << thermostat.p_eta;
        #endif

        std::cout << std::endl; // Si va a capo indipendentemente dal valore degli if
      }

      // Stampa su file
      if (t > pars.eqstps && pars.savemeasure > 0 && (t % pars.savemeasure == 0)){
        save_measures(t);
      }
    } // Fine for sui passi di integrazione

    // Dopo aver finito tutti i passi di integrazione si chiamano i metodi per calcolare e salvare su file la RDF o lo SSF 
    // Ovviamente solo se pars.RDF e/o pars.SSF sono >=1
    if(pars.RDF){
      rad.evaluate_g(pars.Np, pars.rho);
      rad.save_distribution();
    }
    if(pars.SSF){
      ssf.evaluate_S_k(pars.Np);
      ssf.save_ssf();
    }
    

  } // Fine funzione run 
};

#endif