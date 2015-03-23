#include"allheaders.h"
#include<time.h>



// INPUTS
//		step: the maximum distance two successive positions can vary
//		r_init: 2x3 matrix corresponding to the initial positions of the electrons
//		N: the number of positions to be generated

// OUPUT
//		out: a Nx2x3 array of positions

// To use this code, just replace "dummy2" by the function that calculates the wave function

void particle::generate_metropolis(double step)
{
srand(time(NULL));

FILE *pf;
double wave_now, wave_next, condition, rand_num,range;
double x_old;
range=(0.5*RAND_MAX);

pf = fopen("metropol.dat", "w");

for (int i=0 ; i < N ; i++)
{

wave_now=wavefunction();
condition = 0;
rand_num = 0;

rand_num= (double)rand()/(double)RAND_MAX;

x_old=x;




	do{
		x= x_old+ step* (((double)rand()/range)-1); 
		wave_next=wavefunction();
		condition=(wave_next*wave_next)/(wave_now*wave_now) ;

	} while ( condition < rand_num);

fprintf(pf,"%f \n",x);
}

fclose(pf);
}




