#include "general.h"
#include "inicializar.h"
#include "avanzar.h"
#include "Montecarlo.h"
#include <stdio.h>
#include <math.h>

/*tengo dos formas de hacer montecarlo:

1) 
	Muevo las posiciones y momentos de 1 sola partícula y calculo la energía de
	todo el sistema. 
	La comparo con la energía original. 
	Si la energía final es menor que la energía inicial, mantengo el cambio
	Si la energía final es mayor que la energía inicial, hago el cambio si 
	P = e^Delta  E > numero random uniforme.
	
	Ventajas: No sé
	Desventajas: Por cada cambio de 1 sola particula tengo que calcular toda la energía
	(para 1000 partículas son 149 millones de iteraciones)

2)
	Muevo las posiciones y momentos de una sola partícula y calculo la energía 
	para esa sola partícula.
	La comparo con el promedio de energía por particula del sistema (o sea,
	ya calculé la energía de todo el sistema y la dividí por N el número de partículas,
	eso es la energía media por partícula).
	Si la energía final es menor que la energía inicial, mantengo el cambio,
	Si la energía final es mayor que la energía inicial, hago el cambio si la
	P = e^Delta E > numero random uniforme
	
	Ventajas: Menos iteraciones
	Desventajas: no me sirve la implementación que ya tengo que calcula la energía
	
	voy por opcion 1) que está implementado y vemos si se vuelve muy insostenible.
*/

double montecarlo(double* x, double* v, double* X, double *V, double* p_cuad, double eta, double m, double L, int N)
{
	double E_inicial = Energia_paper(x,v,p_cuad,eta,m,L,N);
	double E_final, Delta_E;
	double P_aceptacion, rand;
	double T, T_f, p_f2, p_f;
	int i,j;
	//defino la temperatura de fermi como la energía de fermi (kT_f = e_f, con k = 1)
	//y la energía de fermi a partir del impulso de fermi e_f = p_f^2/2m, entonces
	p_f2 = pow((6*N*pow(M_PI,2))/(L*L*L),2/3.0);
	p_f = pow(p_f2,1/2);
	T_f = p_f2/(2*m);
	//y ahora defino la temperatura T a partir de tau = 0.05 (valor del paper)
	T = 0.05*T_f;
	//escribo un for de 0 N, para todas las partículas
	for(i = 0; i < N; i++){
		//escribo un for de 0 a 3, para las 3 coordenadas de cada particula
		for(j = 0; j < 3; j++){
			//modifico las 3 coordenadas de la partícula en el vector copia (copia(x) = X)(mayuscula)
			X[3*i+j] = random()*L;
			V[3*i+j] = 2*(random()-0.5)*p_f;
			//random()-0.5 me da numeros entre -0.5 y 0.5
			//multiplico por 2 para que me de entre -1 y 1
			//multiplico por p_f para que me de del orden de entre -p_f y p_f
		}
		//Con el vector copia (modificado) calculo la energía del sistema 
		//y calculo el delta E
		E_final = Energia_paper(X,V,p_cuad,eta,m,L,N);
		Delta_E = E_final - E_inicial;
		//escribo el algoritmo de montecarlo: 
		if(Delta_E > 0){
			//si la energía es mayor a 0, defino una probabilidad de aceptacion
			P_aceptacion = exp(-Delta_E/T);
			//defino un número random uniforme para comparar
			rand = random();
			//y si la probabilidad de aceptación es mayor que el número random,
			//entonces, acepto el estado. Esto es, guardar las coordenadas
			//modificadas en el vector original
			if(P_aceptacion > rand){
				for(j = 0; j < 3; j++){
					x[3*i+j] = X[3*i+j];
					v[3*i+j] = V[3*i+j];
				}
			//además, redefino el estado de energía inicial como el estado final
			//porque se aceptó el nuevo estado
			E_inicial = E_final;
			}
			else{
				//si P_aceptacion < rand no se acepta el estado y
				//el vector modificado lo vuelvo a su forma inicial
				E_final = E_inicial;
				for(j = 0; j < 3; j++){
					X[3*i+j] = x[3*i+j];
					V[3*i+j] = v[3*i+j];
				}
			}
		}
		else{
			//Si Delta_E < 0 entonces se acepta directamente.
			//o sea, uso el vector copia modificado para reescribir las coordenadas
			//modificadas del vector original
			//Si no acepto el estado, el vector modificado lo vuelvo a su 
			//estado inicial, copia del vector original x
			E_final = E_inicial;
			for(j = 0; j < 3; j++){
				x[3*i+j] = X[3*i+j];
				v[3*i+j] = V[3*i+j];
			}
		}
	}
	return E_inicial;
}

double montecarlo2(double* x, double* v, double* X, double* V, double* p_cuad, double eta, double m, double L, int N)
{
	double E_inicial = Energia_paper(x,v,p_cuad,eta,m,L,N);
	double E_final, Delta_E;
	double P_aceptacion, rand;
	int i;
	double T, T_f, p_f2;
	//defino la temperatura de fermi como la energía de fermi (kT_f = e_f, con k = 1)
	//y la energía de fermi a partir del impulso de fermi e_f = p_f^2/2m, entonces
	p_f2 = pow((6*N*pow(M_PI,2))/(L*L*L),2/3.0);
	T_f = p_f2/(2*m);
	//y ahora defino la temperatura T a partir de tau = 0.05 (valor del paper)(tau = T/T_f)
	T = 0.05*T_f;
	//escribo un fo
	distribuir_x(X,L,N);
	distribuir_v(V,L,N);
	E_final = Energia_paper(X,V,p_cuad,eta,m,L,N);
	Delta_E = E_final - E_inicial;
	if(Delta_E > 0){
		P_aceptacion = exp(-Delta_E/T);
	}
	else{
		P_aceptacion = 1;
	}
	rand = random();
	if(P_aceptacion > rand){
		for(i = 0; i < 3*N; i++){
			x[i] = X[i];
			v[i] = V[i];
		}
		//además, redefino el estado de energía inicial como el estado final
		//porque se aceptó el nuevo estado
		E_inicial = E_final;
	}
	else{
		//si no se aceptó el cambio, quiero que el vector copia modificado
		//vuelva a ser idéntico al original. Y repito el procedure.
		for(i = 0; i < 3*N; i++){
			X[i] = x[i];
			V[i] = v[i];
		}
		E_final = E_inicial;
	}
	return E_inicial;
}