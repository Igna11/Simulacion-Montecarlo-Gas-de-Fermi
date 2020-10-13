#include "avanzar.h"
#include <stdio.h>
#include <math.h>

double E_cinetica(double* v, double* p_cuad, double m, int N)
{	
	//Calcula el valor de p^2 de cada una de las N partículas
	//Para eso lo que hace es reducir el vector de velocidades de longitud 3N
	//calculando el modulo cuadrado de la velocidad de cada partícula y después
	//agrupándolos en un vector de longitud N llamado p_cuad (malloc)
	int i, k;
	double E_cin = 0;
	//reinicio el vector p_cuad
	for(i = 0; i < N; i++){
		p_cuad[i] = 0;
	}
	for(i = 0; i < N; i ++){
		for(k = 0; k < 3; k++){
			p_cuad[i] += v[3*i+k]*v[3*i+k];
		}
		//voy calculando la energía cinetica
		E_cin += p_cuad[i]/(2*m);
	}
	return E_cin;
}

double Theta(double* p_cuad, double eta, double L, int N)
{
	//Esta función calcula el valor del término del potencial de Pauli del paper
	//donde está la pseudofunción escalón Theta. Para eso lo que hace calcular 
	//la raiz cuadrada de p^2 (sqrt(p_cuad[i])) y con un for voy sumando los valores
	//de la función a una variable Theta
	int i;
	double theta, V_theta, p, p_f, q;
	theta = 0;
	p_f = pow(6*N*pow(M_PI,2)/(L*L*L),1/3.0);
	V_theta = 3.560; //sacado del paper original, pagina 2, tabla 1
	for(i = 0; i < N; i++)
	{
		p = sqrt(p_cuad[i]);
		q = p/p_f;
		theta += V_theta/(1 + exp(-eta*(pow(q,2)-1)));
	}
	return theta;
}
	
double E_potencial_paper(double* x, double* v, double L, int N)	
{	//esta función calcula los términos del potencial de pauli que tienen que ver con
	//pares de partículas. Calcula las distancias entre pares de partículas y la diferencia
	//de momentos entre pares de partículas. Después los mete en una exponencial y los suma.
	int i,j,k;
	double Vp, Vq, r_ij, p_ij, E_pot_paper, r0, p0, p_f;
	E_pot_paper = 0;
	//defino los parámetros que se usan en el paper para reproducir los resultados
	p_f = pow((6*N*pow(M_PI,2))/(L*L*L),1/3.0);
	Vq = 13.517;//sacado del paper original, pagina 2, tabla 1
	Vp = 1.260; //sacado del paper original, pagina 2, tabla 1
	r0 = 0.845/p_f; //sacado del paper más completo, pagina 7, tabla 1
	p0 = 0.193*p_f; //sacado del paper más completo, pagina 7, tabla 1
	//Calcula la suma de las exponenciales e^(r_{ij})
	for(i = 0; i < N; i++)
	{
		for(j = i + 1; j < N; j++)
		{
			for(k = 0; k < 3; k++)
			{
				//Calculo el módulo |x_i - x_(j =/=i)|
				r_ij += abs(x[3*i+k]-x[3*j+k]);
				//calculo el módulo |p_i - p_(j =/=i)|
				p_ij += abs(v[3*i+k]-v[3*j+k]);
			}
			E_pot_paper += Vq*exp(-r_ij/r0) + Vp*exp(-p_ij/p0);
			r_ij = p_ij = 0;
		}
	}
	return E_pot_paper;
}

double Energia_paper(double* x, double* v, double* p_cuad, double eta, double m, double L, int N)
{
	double E_paper = 0;
	double E_cin = 0;
	double E_pot = 0;
	double E_theta = 0;
	E_pot = E_potencial_paper(x,v,L,N);
	E_cin = E_cinetica(v, p_cuad, m, N);
	E_theta = Theta(p_cuad, eta, L, N);
	E_paper = E_pot + E_cin + E_theta;
	return E_paper/N;
}