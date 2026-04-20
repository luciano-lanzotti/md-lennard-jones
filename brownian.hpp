/**
 * @file brownian.hpp
 * @brief Brownian (Langevin) dynamics propagator.
 *
 * Implements the stochastic propagator exp(iL_G) for one particle at a time:
 *
 *   v_new = v * exp(-gamma*dt) + sqrt(T*(1-exp(-2*gamma*dt))/m) * G
 *
 * where G is a 3-component Gaussian vector generated via the Box-Muller transform.
 * The friction coefficient gamma defaults to 1. Temperature is taken from simpars.T_eq.
 */
#ifndef BROWNIAN_HPP
#define BROWNIAN_HPP

#include "./pvector.hpp"
#include "./params.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include <cmath>
#define M_E 2.718281828459
#define M_PI 3.14159265358979323846

template<typename particle_type> class brownian{
public:
    simpars pars;
    double T, gamma, K2, K; // K non è l'energia cinetica
    pvector<double,3> G; // Terna di numeri distribuiti secondo una gaussiana standardizzata

    brownian(){
        T = pars.T_eq; // Temperatura di equilibrio
        gamma = 1; // coefficiente di attrito viscoso
    }

    // Non passo un vettore di particelle, ma la i-esima particella
    // Ad ogni chiamata di questo metodo evolgo 1 particella
    void expiLG(double dt, particle_type& p){
        K2 = T*(1-exp(-2*gamma*dt))*p.m;
        K = sqrt(K2);
        G = {gauss(),gauss(),gauss()};
        p.v = p.v*exp(-gamma*dt) + (K*G/p.m); 
    } 

    // Generatore di numeri gaussiani secondo box muller
    double gauss(void){
        double x1, x2;
        do{
            x1 = rng.ranf();
            x2 = rng.ranf();
        }while (x1==0 || x2==0);
        return cos(2*M_PI*x2)*sqrt(-2.0*log(x1));
    }

   


    
};

#endif
