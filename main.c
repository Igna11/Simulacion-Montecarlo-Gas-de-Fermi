/* gcc -Wall -O3 -o FermiGas.exe main.c general.c inicializar.c avanzar.c Montecarlo.c -lm*/
#include "general.h"
#include "inicializar.h"
#include "avanzar.h"
#include "Montecarlo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int main(int argc, char *argv[]){
	// Definicion del main
	int N, iteraciones, i;
	double rho, rho_i, rho_f, paso, L, eta, m;
	eta = 30;
	m = 1;

//----------------------CONTADOR--------------------------------------//
	time_t timer;
    char buffer[26];
    struct tm* tm_info;

    timer = time(NULL);
    tm_info = localtime(&timer);
    

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
	double *X = (double *) malloc(3*N*sizeof(double)); //aca guardo los valores modificados de las coordenadas
	double *V = (double *) malloc(3*N*sizeof(double)); //aca guardo los valores modificados de las coordenadas

	
//-------REINICIAR LOS MALLOCS ES UNA BUENA PRÁCTICA ESTANDAR---------//
	for(rho = rho_i; rho < rho_f; rho = rho + paso){

	for(i = 0; i < 3*N; i++){
		x[i] = 0;
		v[i] = 0;
		X[i] = 0;
		V[i] = 0;
 	}
	for(i = 0; i < N; i++){
		p_cuad[i] = 0;
	}

	
//-------------------DISTRIBUYO LAS PARTÍCULAS------------------------//

	
	//defino el cubo donde trabajo
	//defino el cubo donde trabajo
	// L = raizcubica(N/rho).
	L = cbrt(N/rho);
	
	//Corro las funciones que fui armando
	distribuir_x(x,L,N); // distribuir_x está en inicializar.c
	distribuir_v(v,L,N); // distribuir_v está en inicializar.c
	
	//genero una copia de los vectores originales de posiciones y momentos
	for(i = 0; i < 3*N; i++)
	{
		X[i] = x[i];
		V[i] = v[i];
	}
	

//---------------------NOMBRE DEL TXT SEGUN PARÁMETROS----------------//
	char filename[255];
//--------------------- ARRANCO CON EL PROGRAMA ----------------------//


	//inicializo la variable para la energía
	FILE* fp;
	
	double Ener;
	double E_cin;
	
		sprintf(filename, "02 - Ecin_vs_rho - N%i - rho %f.txt",N, rho);
		fp = fopen(filename, "w");
		for(i = 0; i < iteraciones; i++){
			Ener = montecarlo(x,v,X,V,p_cuad,eta,m,L,N);
			E_cin = E_cinetica(v, p_cuad, m, N);
			printf("rho = %f, iteracion N%i, Energia cinetica = %f Energia = %f\n", rho, i, E_cin ,Ener);
			fprintf(fp,"%f\t%f\n", E_cin, Ener);
		}
		fclose(fp);
	}
	
	//----------------------TIMER --------------------------//
	
	
	strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
    puts(buffer);//Imprime la hora a la que se empezó a correr la simulación
	timer = time(NULL);
    tm_info = localtime(&timer);
	strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
    puts(buffer);//Imprime la hora a la que termino de correr la simulación
	
	
	//------------------------------------------------------------//
	//printf("energia paper configuracion N%i= %f",iteraciones, E_paper); 
	
	/*//-- LINEA PARA EXPORTAR .txt QUE SIEMPRE MEOVLIDO COMO SE ESCRIBE --
	FILE* fp;
	fp = fopen("cinetica.tx
	for(i = 0; i < N; i++)
	{
		fprintf(fp,"%f\t",p_cuad[i]);
		for(k = 0; k < 3; k++)
		{
			fprintf(fp, "%f\t", v[3*i+k]);
		}
		fprintf(fp, "\n");
	}
	*/
	//libero el espacio de memoria
	free(x);
	free(v);
	free(p_cuad);
	free(X);
	free(V);
	return 0;
}