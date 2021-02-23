/* 
 * File:   Simulation.h
 * Author: Charlie
 *
 * Created on 08 February 2016, 12:35
 */

#ifndef SIMULATION_H
#define	SIMULATION_H

#include "General.h"
#include "Analysis.h"

void Run(int pot, double m, double c, int nP, double dx, int nS, double dt, int nPrint, int nAnalysis);

pair< Sim, vector< complex<double> > > initialise(int nP, int nS, double dr, double dt, double m, int Pot, double C, string name, int nPrint, int nAnalysis);

vector< complex<double> > appM1M2(vector< complex<double> > U, Sim sim);
vector< complex<double> > appM2M1(vector< complex<double> > U, Sim sim);

vector<double> V(vector< complex<double> > U, Sim sim);
vector<double> fastSNV(vector< complex<double> > U, Sim sim);
vector< complex<double> > appExpVU(vector< complex<double> > U, Sim sim);

vector< complex<double> > Update(vector< complex<double> > U, Sim sim);

vector< complex<double> > exact(Sim sim, double t);

void printSimDetails(Sim sim, string name);

void setupPrintS(string name, int nP, double dx);
void printSim(vector< complex<double> > U, string name, double t);

#endif	/* SIMULATION_H */

