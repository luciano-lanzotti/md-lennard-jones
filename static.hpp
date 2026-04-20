/**
 * @file static.hpp
 * @brief Static structure factor S(k) accumulator and evaluator.
 *
 * S(k) = (1/N) * <|rho(k)|^2>  where  rho(k) = sum_j exp(i k · r_j)
 *
 * At each production step update_ssf() randomly draws `total_points` wave vectors
 * k = (2pi/L) * n  with n a random integer triple in [-n_max, n_max]^3,
 * accumulates |rho(k)|^2 in a histogram binned by |k|, and tracks the hit count.
 *
 * Usage pattern:
 *   1. Call init_ssf(L) after the box size is known.
 *   2. Call update_ssf(parts, Np, L) at every production step.
 *   3. Call evaluate_S_k(Np) and save_ssf() at the end of the run.
 */
#ifndef STATIC_HPP
#define STATIC_HPP

#include "./pvector.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include <cmath>
#include <string>
#include <complex>
#include <cstdlib>

// Creo una classe per calcolare lo static structure factor
template<typename particle_type> class ssf{
    long int n_bins; // Numero di bin del dominio di k
    long int total_points; // Numero di vettori che generiamo ad ogni chiamata di update_ssf
    double bin_size; // Dimensione di un bin. In generale diversio da dk
    double k_max; // Massimo valore del modulo del vettore k
    double n_max; // Massimo intero nella formula del vettore k
    std::vector<double> count; // Array che conterrà quante volte abbiamo sorteggiato un modulo k
    std::vector<double> rho2; // Modulo quadro di rho
    std::vector<double> S_k; // statick structure factor
    std::string output_name; // Nome del file dove salveremo k e S_k
    public:
    ssf(){
        total_points = 1000; 
        n_max=15;
        output_name = "dati/ssf/ssf.dat";    
    }

    void init_ssf(pvector<double,3> L){
        // I vettori con il modulo, non nullo, più piccolo possibile sono: (2pi/L, 0, 0); (0, 2pi/L, 0); (0, 0, 2pi/L)
        // Allora il modulo più piccolo è dk = 2pi/L. Tuttavia non possiamo prendere dk come stima del binsize.
        // Infatti supponiamo i vettori (dk, dk, dk) e (dk, dk, 0) i loro moduli sono rispettivamente dk*sqrt(3) e dk*sqrt(2) che si differenziano per meno di dk.
        // In conclusione se scegliessimo come binsize dk, staremo sovrastimando la dimensione del bin. 
        // Per convenzione allora binsize è una frazione di dk.  
        double dk = 2*M_PI/((L(0)+L(1)+L(2))/3.0); // Facciamo una media di L, tanto dk serve solo come stima.
        bin_size = dk/3.0;
        pvector<double,3> k_max_vec = {1/L(0), 1/L(1), 1/L(2)}; // Costruiamo un vettore che nella i-esima componente ha 1/L_i
        k_max_vec *= 2*M_PI*n_max; // Questo è esattamente il più 'grande' vettore che possiamo costruire dati n_max 
        k_max = k_max_vec.norm(); // Ne calcoliamo la norma
        n_bins = (long int)rint(k_max/bin_size)+1; // Calcoliamo quanti bin servono. Si ricorda che il binsize è funzione di dk. 
        // Gli array rho2, S_k e count hanno dimensione n_bins e li inizializziamo a 0.
        rho2.resize(n_bins);
        S_k.resize(n_bins);
        count.resize(n_bins);
        for(int i=0; i<n_bins; i++){
            rho2[i]=0.;
            S_k[i]=0.;
            count[i]=0;
        }
    }

    
    void update_ssf(std::vector<particle_type>& parts, int Np, pvector<double,3> L){
        std::complex<double> rho; // rho(k) definito come sum_j exp(i*kj*rj). Sarà un numero complesso.
        pvector<double,3> n; // Vettore che conterrà la terna di n, per convenienza creaiamo di tipo double. 
        pvector<double,3> k_vec; // pvector che conterrà il vettore k generato casualemente.
        double k; // Modulo del vettore casuale.
        long int bin; // Indice del bin a cui k appartiene
        for(int p=0; p<total_points; p++){ // For sui 2000 vettori casuali
            // Generiamo una terna casuale di numeri interi nx,ny,nz.
            for(int i=0; i<3; i++){
                n(i) = (double)random_int(-n_max,n_max);
            }
            k_vec = n.divcw(L)*2*M_PI; // Dalla terna di n otteniamo un vettore
            k = k_vec.norm(); // Ne calcoliamo la norma
            bin = (long int)rint(k/bin_size); // Calcoliamo a quale bin corrisponde quella norma
            rho = evaluate_rho(k_vec, parts, Np); // Calcoliamo rho associato a quel vettore
            rho2[bin] += std::norm(rho); // Aggiorniamo rho2 nella posizione bin-esima
            count[bin] += 1; // Aggiorniamo coount nella posizione bin-esima
        }    
    }

    // Funzione per generare un interno nell'intevallo [min,max]. Entrambi gli estremi inclusi
    long int random_int(long int min, long int max){
        long int r = max-min+1;
        return rand()%r + min;
    }
    
    // Funzione per calcolare rho(k)
    std::complex<double> evaluate_rho(pvector<double,3> k_vec, std::vector<particle_type>& parts, int Np){
        std::complex<double> rho;
        double arg;
        // Inizializziamo rho(k) a zero. Prima di iniziare a sommare su tutti i possibili contributi e^(i*kj*rj)  
        rho = {0.0, 0.0}; 
        for(int i=0; i<Np; i++){
            arg = k_vec * parts[i].r;
            rho += std::polar(1.,arg);  // Creaiamo un numero complesso del tipo e^{i*kj*rj}. Con j=1,...,N
        }
        return rho;
    } 

    // Calcoliamo il fattore di struttura
    void evaluate_S_k(long int Np){
        for(int bin=0; bin<n_bins; bin++){
            // Se il bin è vuoto sono possibili due casi:
            // 1. Poiché i vettori k sono discreti, non è stato possibile generare un vettore con modulo in quel range.
            // Ad esempio è impossibile generare un vettore non nullo con modulo minore di dk=2pi/L.
            // 2. Per caso non sono stati estratti vettori con modulo in quel range. Tipicamente ciò succede per k molto grandi. 
            
            // Ponendo S_k[bin] a 0 se count[bin] è nullo, evitiamo una divisione tra due 0.
            if (count[bin]==0){
                S_k[bin]=0.0;
            }
            // rho2[bin]/count[bin] = <rho2[bin]>. Allora dividendo anche per il numero di particelle, si ottiene S_k
            else{
                S_k[bin] = rho2[bin]/(Np*(double)count[bin]);
            } 
            
        }
    }

    void save_ssf(){
        std::fstream f;
        f.open(output_name,std::ios::out|std::ios::trunc); 
        for(int n=0; n<n_bins; n++){
            f << bin_size*(n+0.5) << "\t ";  // Salviamo il valore centrale del bin
            f << S_k[n] << "\n ";
        }
        f.sync();      
        f.close();
    }
    // Funzione che permette di cambiare nome al file di output. Utile per non sovrascrivere files.
    void set_output_name(std::string new_name){
        output_name = new_name;
    }


};
#endif