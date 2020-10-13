#ifndef MONTECARLO_H
#define MONTECARLO_H
#include "math.h"

double montecarlo(double* x, double* v, double* X, double *V, double* p_cuad, double eta, double m, double L, int N);
double montecarlo2(double* x, double* v, double* X, double *V, double* p_cuad, double eta, double m, double L, int N);

#endif