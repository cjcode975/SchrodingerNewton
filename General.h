/* 
 * File:   General.h
 * Author: Charlie
 *
 * Created on 08 February 2016, 12:35
 */

#ifndef GENERAL_H
#define	GENERAL_H

#include <cstdlib>
#include <complex>
#include <cmath>
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>

using namespace std;

struct Sim {    
    int nPnts, nStps, pot, nPrint, nAnalysis;
    double M, dr, dt, W;
    complex<double> aa, ab, bb, ac, bc;
};

double normSq(complex<double> c);
complex<double> bar(complex<double> c);

#endif	/* GENERAL_H */

