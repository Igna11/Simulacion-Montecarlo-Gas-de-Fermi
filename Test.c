#include "avanzar.h"
#include "Montecarlo.h"
#include <stdio.h>
#include <stdlib.h>


double test_distribuir(double* x, 
						double* v, 
						double* X, 
						double* V, 
						double* p_cuad,
						double eta,
						double m, 
						double L, 
						int N,
						int l){
	/*
		Esta función se encarga de imprmir en .txt las posiciones de las partículsa
		para ver si se está respetando la condición periódica de contorno
	*/
	int i,j;
	FILE* fp;
	char filename[255];
	sprintf(filename, "C:/Users/igna/Desktop/Igna/Facultad/Fisica computacional/FINAL/Codigo/src/Tests/CCP/TEST_CCP_%i.txt",l);
	fp = fopen(filename, "w");
	for(i = 0; i < N; i++){
		for(j = 0; j < 3; j++){
			fprintf(fp,"%f\t", x[3*i + j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return 0;
}