/* 
 * File:   Analysis.h
 * Author: Charlie
 *
 * Created on 08 February 2016, 10:59
 */

#ifndef ANALYSIS_H
#define	ANALYSIS_H

#include "General.h"
#include "Simulation.h"

pair< vector<double>, vector< complex<double> > > analysis(vector< complex<double> > U, Sim sim, double t);

void setupPrintA(int nAn, int nS, double dt);
void printAnalysis(vector< complex<double> > U, Sim sim, double t, double m, bool m_first);

#endif	/* ANALYSIS_H */

