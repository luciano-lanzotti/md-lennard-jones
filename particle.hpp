/**
 * @file particle.hpp
 * @brief Particle base class and Lennard-Jones particle derived class.
 *
 * `particle` stores position r, velocity v, force f, and LJ parameters.
 * Provides symplectic propagators expiLp (momentum half-step) and expiLq (position step).
 *
 * `particleLJ` adds:
 *   - vij()  : LJ pair potential (no shift)
 *   - vijs() : LJ pair potential shifted by vcut
 *   - fij()  : LJ pair force, respecting periodic boundary conditions
 *   - set_vcut() : precomputes the potential value at rc for shifted interactions
 *
 * Both classes apply minimum-image periodic boundary conditions via the PBC formula
 *   Dr -= L * rint(Dr / L)
 */
#ifndef _PARTICLE_
#define _PARTICLE_

#include "./pvector.hpp"
#include "./params.hpp"
using ntype = double;

class particle{
protected:
  ntype vcut;
public:
  ntype sigma, epsilon, rc, m;
  pvector<ntype,3> r, v, f; 
  void set_vcut(void);
  void set_sigma(ntype sig){
    sigma = sig;
  }
  void set_epsilon(ntype eps){
    epsilon = eps;
  }
  void set_rcut(ntype rcut){
    rc = rcut;
  }
  particle(){
    simpars pars; 
    sigma = pars.sigma;
    epsilon = pars.epsilon;
    rc = pars.rc;
    m = pars.mass;
    vcut = 0.0;
  }

  // Evolutori per dinamica molecolare
  void expiLp(ntype dt){
    v += (f/m)*dt;
  }
  void expiLq(ntype dt){
    r += v*dt;
  }

  
};

class particleLJ: public particle{
public:
  // Metodo per calcolare l'interazione con cut-off non riscalata
  ntype vij(particleLJ P, pvector<ntype,3> L){
    ntype ene;
    pvector<ntype,3> Dr;
    Dr = r - P.r;
    Dr = Dr - L.mulcw(rint(Dr.divcw(L))); // Usiamo le PBC
    ntype rn = Dr.norm();
    ntype rsq = rn*rn; 
    if (rsq < rc*rc) // Interazione non nulla solo se la distanza è minore del raggio di cut-off
      ene = 4.0*epsilon*(pow(sigma/rn,12.0)-pow(sigma/rn,6));
    else
      ene = 0;
    return ene;
  }

  // Metodo per calcolare l'interazione con cut-off e riscalata
  ntype vijs(particleLJ P, pvector<ntype,3> L){
    ntype ene;
    pvector<ntype,3> Dr;
    Dr = r - P.r;
    Dr = Dr - L.mulcw(rint(Dr.divcw(L))); // Usiamo le PBC
    ntype rn = Dr.norm();
    ntype rsq = rn*rn; 
    if (rsq < rc*rc) // Interazione non nulla solo se la distanza è minore del raggio di cut-off
      ene = 4.0*epsilon*(pow(sigma/rn,12.0)-pow(sigma/rn,6)) - vcut;
    else
      ene = 0;
    return ene;
  }

  // Calcolo del potenziale al cut-off
  void set_vcut(void){
    ntype rijsq, srij2, srij6, srij12;
    ntype epsilon4=epsilon*4.0;
    rijsq = rc*rc;
    srij2 = sigma*sigma / rijsq;
    srij6  = srij2 * srij2 * srij2;
    srij12 = srij6*srij6;
    vcut = epsilon4*(srij12 - srij6);
  }

  // Forza che agisce tra la particella di questa classe e un'altra passata in argomento
  pvector<ntype,3> fij(particleLJ P, pvector<ntype,3> L){
    pvector<ntype,3> fijv;
    pvector<ntype,3> Dr = r-P.r;
    Dr = Dr - L.mulcw(rint(Dr.divcw(L))); 
    ntype rn = Dr.norm();
    ntype rsq = rn*rn;
    if (rsq < rc*rc){
      fijv =  4 * epsilon * ( (12*pow(sigma/rn,12)) - (6*pow(sigma/rn,6))) * (Dr/rsq); // Forza ottenuta con la meno derivata del potenziale di LJ
    }
    else fijv = {0.,0.,0.}; // Se le particelle sono più lontane di rc, allora la forza tra di loro è 0
    return fijv;
  }

};

#endif 
