#ifndef AVANZAR_H
#define AVANZAR_H
#include "math.h"

double E_cinetica(double* v, double* p_cuad, double m, int N);
double Theta(double* p_cuad, double eta, double L, int N);
double E_potencial_paper(double* x, double* v, double L, int N);
double Energia_paper(double* x, double* v, double* p_cuad, double eta, double m, double L, int N);
#endif