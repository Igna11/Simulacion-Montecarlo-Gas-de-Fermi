/* gcc -Wall -O3 -o FermiGas.exe main.c general.c inicializar.c avanzar.c Montecarlo.c -lm*/
#include "general.h"
#include "inicializar.h"
#include "avanzar.h"
#include "Montecarlo.h"
#include "Test.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int main(int argc, char *argv[]){
	// Definicion del main
	int N, iteraciones, i;
	double rho, rho_i, rho_f, paso, L, eta, m, p_f, p_f2, E_f, timer1, timer2, t1, t2;
	eta = 30;
	m = 1;

//----------------------CONTADOR--------------------------------------//

    timer1 = time(NULL);

//---------------------- PARAMETROS INPUT ----------------------------//


	printf("\nPasame el numero de particulas ameo\n");
	scanf("%int", &N);
	printf("\nPasame el numero de iteraciones\n");
	scanf("%i", &iteraciones);
	printf("\nPasame la densidad inicial rho_i\n");
	scanf("%lf", &rho_i);
	printf("\nPasame la densidad final rho_f\n");
	scanf("%lf", &rho_f);
	printf("\nDame el paso entre densidades\n");
	scanf("%lf", &paso);


	//srand(time(0)) pone una semilla inicial que varía con el tiempo para 
	//maximizar la irreproducibilidad del conjunto de numeros randoms generados
	srand(time(0));

	//Tengo el numero de partículas N, tengo la temperatura T, tengo la densidad rho
	//Ahora quiero distribuir estas partículas de forma random en un recinto de lado



//------------------------ESPACIOES DE MEMORIA------------------------//


	//defino los espacios de memoria para vectores: 3*N para 3 coordenadas por partícula
	double *x = (double *) malloc(3*N*sizeof(double)); //posición de cada partícula
	double *v = (double *) malloc(3*N*sizeof(double)); //velocidad de cada partícula
	double *p_cuad = (double *) malloc(N*sizeof(double)); //modulo cuadrado del momento de cada partícula
	double *x_modif = (double *) malloc(3*N*sizeof(double)); //aca guardo los valores modificados de las coordenadas
	double *v_modif = (double *) malloc(3*N*sizeof(double)); //aca guardo los valores modificados de las coordenadas


//-------REINICIAR LOS MALLOCS ES UNA BUENA PRÁCTICA ESTANDAR---------//
	for(rho = rho_i; rho < rho_f; rho = rho + paso){

		for(i = 0; i < 3*N; i++){
			x[i] = 0;
			v[i] = 0;
			x_modif[i] = 0;
			v_modif[i] = 0;
		}
		for(i = 0; i < N; i++){
			p_cuad[i] = 0;
		}


	//-------------------DISTRIBUYO LAS PARTÍCULAS------------------------//


		//defino el cubo donde trabajo
		//defino el cubo donde trabajo
		// L = raizcubica(N/rho).
		L = cbrt(N/rho);
		//calculo el impulso de fermi p_f, el impulso de fermi cuadrado p_f2
		//la energía de fermi E_f = p_f2/(2*m) asi no lo tiene que hacer cada vez
		//que calculo una energia. Veremos cuanto mejora el rendimiento (01/11/2020)
		p_f = pow(6*pow(M_PI, 2)*rho, 1/3.0);
		p_f2 = pow(p_f, 2.0);
		E_f = p_f2/(2*m);

		//Corro las funciones que fui armando
		distribuir_x(x,L,N); // distribuir_x está en inicializar.c
		distribuir_v(v,p_f,L,N); // distribuir_v está en inicializar.c

		//genero una copia de los vectores originales de posiciones y momentos
		for(i = 0; i < 3*N; i++)
		{
			x_modif[i] = x[i];
			v_modif[i] = v[i];
		}
		

	//---------------------NOMBRE DEL TXT SEGUN PARÁMETROS----------------//
		char filename[255];
	//--------------------- ARRANCO CON EL PROGRAMA ----------------------//


		//inicializo la variable para la energía
		FILE* fp;

		double Ener;
		double p, q;

		sprintf(filename, "05 - Ecin_vs_rho - N%i - rho %f.txt",N, rho);
		fp = fopen(filename, "w");
		t1 = time(NULL);
		for(i = 0; i < iteraciones; i++){
			Ener = montecarlo(x,v,x_modif,v_modif,p_cuad,eta,m,L,rho,p_f,p_f2,E_f,N);
			p = pow(2*m*E_cinetica(v, p_cuad, m, N), 1/2.0);
			q = p/p_f;
			//printf("rho = %f, iteracion N%i, q = p/p_f = %f Energia = %f\n", rho, i, q ,Ener);
			fprintf(fp,"%f\t%f\n", q, Ener);
		}
		t2 = time(NULL) - t1;
		printf("\ntiempo por cada densidad: rho = %f, t = %f\n", rho, t2);
		fclose(fp);
	}

	//----------------------TIMER --------------------------//

	timer2 = time(NULL) - timer1;
	printf("\nTiempo total: %f\n", timer2);

	//------------------------------------------------------------//

	//libero el espacio de memoria
	free(x);
	free(v);
	free(p_cuad);
	free(x_modif);
	free(v_modif);
	return 0;
}