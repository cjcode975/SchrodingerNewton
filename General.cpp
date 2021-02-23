#include "General.h"

double normSq(complex<double> c){ return c.real()*c.real() +c.imag()*c.imag(); }
complex<double> bar(complex<double> c){ return complex<double>(c.real(),-1*c.imag()); }