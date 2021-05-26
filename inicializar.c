#include "general.h"
#include "inicializar.h"
#include "stdlib.h"
#include "math.h"

double distribuir_x(double *x, double L, int N)
{
	//distribuye las partículas en un espacio de dimensiones LxLxL
	//con L = raizcubica(N/rho), generando numeros random que van
	//de 0 a L y se acomodan en un espacio de memoria de dimensión
	//3N
	int i;
	for(i = 0; i < 3*N; i++)
	{
		x[i] = random()*L;
	}
	return 0;
}
double distribuir_v(double *v, double p_f, double L, int N)
{
	/*distribuye las componentes de las velocidades de forma uniforme
	solo para la inicialización del sistema. Hay que darle la masa
	Tener en cuenta que van a ser todas velocidades positivas*/
	int i;

	for(i = 0; i < 3*N; i++)
	{
		//Distribuyo las velocidades con un valor random uniforme entre
		//0 y el impulso de Fermi.
		//random()-0.5 me da numeros entre -0.5 y 0.5
		//multiplico por 2 para que me de entre -1 y 1
		//multiplico por p_f para que me de del orden de entre -p_f y 
 		v[i] = 2*(random()-0.5)*p_f;
	}
	return 0;
}
