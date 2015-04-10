#include "allheaders.h"
#include <stdio.h>

using namespace std;

int main(void)

{

////////////////  PARAMETERS   /////////////////

double metro_step = 2;
double s_min = 0.4;
double s_max = 5;
double s_step = 0.02;
int it = 10000; 		// Number of iterations in the Monte Carlo 
int N_beta = 10;

////////////////   ALGORITHM   /////////////////

int N = (int)((s_max - s_min)/s_step); 	// Number of values of s
double pos[6]= {1,1,1,1,1,1};	// Initial position of electrons
double* min; 			
particle electron(pos,2.1,0.5,it);	// Don't care about the arguments, just to initialize
electron.wavefunction(); 		// Initializes the variables in electron
double s,energy,beta;
FILE* data = fopen("test.dat","w");

for (int j=0 ; j <= N+1; j++)
{	
	s=s_min+s_step*j;
	energy = 0;beta = 0;
	for (int i = 0; i < 10; ++i)
	{
		if (i>4)
		{
			min = minimize(electron,s, N_beta);
			energy += min[0]/5;
		}
		
	}
	 
	fprintf(data, "%f %f\n", s,energy);
	fprintf(stdout,"%2.1f percent\n", 100*(float)j/(float)(N+1));
}

fclose(data);
return 0;
}
