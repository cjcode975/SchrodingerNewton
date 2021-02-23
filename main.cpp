/* 
 * File:   main.cpp
 * Author: Charlie
 *
 * Created on 08 February 2016, 10:58
 */

#include "General.h"
#include "Simulation.h"

int main(int argc, char** argv) {

    int pot = 2; //potential type- 0 for free, 1 for harmonic, 2 for SN
    double w = 0.25;
    
    double start_mass = 1.0; 
    double inc_mass = 0.1;
    int nInc = 1;
    
    int nP = 1000;
    double dx = 0.25;
    
    int nS = 600;
    double dt = 0.0005;
    
    int nPrint = 50;
    int AeNS = 200;
    
    setupPrintA(AeNS, nS, dt);
    
    for (int i=0; i<nInc; i++){
        double mass = start_mass+i*inc_mass;
        Run(pot, mass, w, nP, dx, nS, dt, nPrint, AeNS);
        cout << "Done mass " << mass << endl;
    }        
    
    return 0;
}

