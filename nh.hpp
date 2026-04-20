/**
 * @file nh.hpp
 * @brief Nosé-Hoover chain thermostat.
 *
 * Implements the extended-system Nosé-Hoover equations of motion split into
 * three Liouville propagators suitable for a Trotter-factored integrator:
 *
 *   expiLeta(dt)         -- update thermostat coordinate eta
 *   expiLpeta(dt, parts) -- update thermostat momentum p_eta from kinetic energy
 *   expiLFNH(dt, parts)  -- rescale all particle velocities by exp(-p_eta * dt / Q)
 *
 * The thermal mass is Q = Nf * T * tau^2 where tau is the LJ characteristic time.
 * energy() returns the conserved thermostat contribution P_eta^2/(2Q) + Nf*T*eta.
 */
#ifndef THERMOSTAT_HPP
#define THERMOSTAT_HPP

#include "./pvector.hpp"
#include "./params.hpp"
#include "./particle.hpp"
#include <cmath>
#define M_E 2.718281828459

template<typename particle_type> class nh {
public:
    double eta;         
    double p_eta;       
    double Q;
    double tau; // tempo tipico del sistema
    double T; // Temperatura del sistema
    simpars pars;   // Parametri della simulazione        
    nh(){
        // inizializziamo a 0 eta e p_eta
        eta=0; 
        p_eta=0;
        tau = pars.sigma*sqrt(pars.mass/pars.epsilon); // Il tempo caratteristico è dettato dai parametri della simulazione
        T = pars.T_eq; // La temperatura T rappresenta quella che vorremmo mantenere, quindi scegliamo la T_eq dalla classe parametri 
        Q = pars.Nf * T * tau * tau; // Con questa formula si calcola Q
    }

    // Evolutore variabile eta
    void expiLeta(double dt){
        eta += (p_eta/Q)*dt;
    }
    
    // Evolutore variabile p_eta
    //Passiamo per riferimento, per evitare di copiare ogni volta il vettore di particelle
    void expiLpeta(double dt, std::vector<particle_type>& parts){
        double F = 0;
        // for su tutte le particelle
        for (int i=0; i<pars.Np; i++){
            F+=parts[i].m * ( parts[i].v * parts[i].v );
        }
        F -= pars.Nf*T; // pars.Nf è il numero di gradi di libertà del sistema e vale 3N
        p_eta += F*dt; // Alla fine p_eta -> p_eta + F*dt
    }

    // Evolutore expiLFNH
    void expiLFNH(double dt, std::vector<particle_type>& parts){
        for (int i=0; i<pars.Np; i++){
            parts[i].v = parts[i].v * exp(-p_eta*dt/Q);
        }
    }
    // Parte dell'energia conserva di nosè hoover che dipende dalle variabili del termostato
    double energy(void){
        return p_eta*p_eta/(2*Q) + (pars.Nf*T*eta);
    }


    
};

#endif
