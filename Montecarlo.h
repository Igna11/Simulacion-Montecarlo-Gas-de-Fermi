#ifndef MONTECARLO_H
#define MONTECARLO_H
#include "math.h"

double montecarlo(double* x,
					double* v,
					double* x_modif,
					double* v_modif,
					double* p_cuad, 
					double eta,
					double m,
					double L,
					double rho,
					double p_f,
					double p_f2,
					double E_f,
					int N);
/*double montecarlo2(double* x,
					double* v,
					double* x_modif,
					double* v_modif,
					double* p_cuad,
					double eta,
					double m,
					double L,
					double rho,
					int N);
*/
#endif