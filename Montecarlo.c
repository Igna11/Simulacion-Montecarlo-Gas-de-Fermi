#include "general.h"
#include "inicializar.h"
#include "avanzar.h"
#include "Montecarlo.h"
#include <stdio.h>
#include <math.h>

double montecarlo(double* x, double* v, double* x_modif, double *v_modif, double* p_cuad, double eta, double m, double L, double rho, double p_f, double p_f2, double E_f, int N)
{
	double E_inicial = Energia_paper(x,v,p_cuad,eta,m,L,rho,p_f,N);
	double E_final, Delta_E;
	double P_aceptacion, rand;
	double T, T_f;	
	int i,j;
	//defino la temperatura de fermi como la energía de fermi (kT_f = e_f, con k = 1)
	//y la energía de fermi a partir del impulso de fermi e_f = p_f^2/2m, entonces
	T_f = E_f;
	double salto_x = L; //defino el salto máximo posible de la partícula como la longitud característica de cada partícula
	double salto_v = p_f/m; //defino el salto máximo posible de la partícula como p_f característica de cada partícula
	//y ahora defino la temperatura T a partir de tau = 0.05 (valor del paper)
	T = 0.05*T_f;
	//escribo un for de 0 N, para todas las partículas
	for(i = 0; i < N; i++){
		//escribo un for de 0 a 3, para las 3 coordenadas de cada particula
		for(j = 0; j < 3; j++){
			//modifico las 3 coordenadas de la partícula en el vector copia 
			//modificación: 18-10-2020 03:05: Ahora voy a "correr" aleatoriamente la partícula, 
			//es decir, dada su posición actual, le sumo a cada coordenada un valor random cuya
			//magnitud voy a tener que explorar con paciencia. Idem para momentos
			x_modif[3*i + j] = x[3*i + j] + (random() - 0.5)*salto_x;
			v_modif[3*i + j] = v[3*i + j] + (random() - 0.5)*salto_v;
			//ahora meto condición periódica de contorno (a las coordenadas espaciales)
			if(x_modif[3*i + j] > L){
				x_modif[3*i + j] = x_modif[3*i + j] - L;
			}
			else if(x_modif[3*i + j] < 0){
				x_modif[3*i + j] = x_modif[3*i + j] + L;
			}
		}
		//Con el vector copia (modificado) calculo la energía del sistema 
		//y calculo el delta E
		E_final = Energia_paper(x_modif,v_modif,p_cuad,eta,m,L,rho,p_f,N);
		Delta_E = E_final - E_inicial;

//----------------------MONTECARLO-----------------------//

		if(Delta_E > 0){
			//si la energía es mayor a 0, defino una probabilidad de aceptacion
			P_aceptacion = exp(-Delta_E/T);
			//defino un número random uniforme para comparar
			rand = random();
			//y si la probabilidad de aceptación es mayor que el número random,
			//entonces, acepto el estado. Esto es, redefinir la energía inicial
			//como la energía final para usarla en el proximo paso y modificar 
			//el vector original para que sea igual al modificado aceptado
			if(P_aceptacion > rand){
				E_inicial = E_final;
				for(j = 0; j < 3; j++){
					x[3*i + j] = x_modif[3*i + j];
					v[3*i + j] = v_modif[3*i + j];
				}
			}
			else{
				//si P_aceptacion < rand no se acepta el estado y
				//el vector modificado lo vuelvo a su forma inicial
				for(j = 0; j < 3; j++){
					x_modif[3*i + j] = x[3*i + j];
					v_modif[3*i + j] = v[3*i + j];
				}
			}
		}
		else{
			//si se acepta el estado, la energia inicial tiene que ser igual a la final, para
			//usarla en el siguiente paso
			E_inicial = E_final;
			for(j = 0; j < 3; j++){
				x[3*i + j] = x_modif[3*i + j];
				v[3*i + j] = v_modif[3*i + j];
			}
		}
	}
	return E_inicial;
}
//////////-------------------
//////////-------------------
//////////-------------------
//////////-------------------


/*double montecarlo2(double* x, double* v, double* x_modif, double* v_modif, double* p_cuad, double eta, double m, double L, double rho, int N)
{
	double E_inicial = Energia_paper(x,v,p_cuad,eta,m,L,rho,N);
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
	distribuir_x(x_modif,L,N);
	distribuir_v(v_modif,L,N);
	E_final = Energia_paper(x_modif,v_modif,p_cuad,eta,m,L,rho,N);
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
			x[i] = x_modif[i];
			v[i] = v_modif[i];
		}
		//además, redefino el estado de energía inicial como el estado final
		//porque se aceptó el nuevo estado
		E_inicial = E_final;
	}
	else{
		//si no se aceptó el cambio, quiero que el vector copia modificado
		//vuelva a ser idéntico al original. Y repito el procedure.
		for(i = 0; i < 3*N; i++){
			x_modif[i] = x[i];
			v_modif[i] = v[i];
		}
		E_final = E_inicial;
	}
	return E_inicial;
}*/