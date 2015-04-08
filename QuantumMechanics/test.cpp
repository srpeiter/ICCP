#include "allheaders.h"
#include <stdio.h>

using namespace std;

int main(void)

{

////////////////  PARAMETERS   /////////////////

double metro_step = 0.1;
double s_min = 0.2;
double s_max = 0.8;
double s_step = 0.01;
int it = 100000; 		// Number of iterations in the Monte Carlo 
double beta_0 = 8.8; 		// Begining values for the minimiztion algorithm
double beta_1 = 9;
double conv_indicator = 0.05;	// The minimization algorithm stops when the difference between two successive values is lower
double eps = 0.1*beta_0; 	// Speed of convergeance in minimization

////////////////   ALGORITHM   /////////////////

int N = (int)((s_max - s_min)/s_step); 	// Number of values of s
double pos[2][3]= {{1,1,1},{1,1,1}};	// Initial position of electrons
double *out; 				// Used to store the minimization results: the minimal energy and the beta at which this is achieved
double dat1[N], dat2[N]; 		// Stores the s values and the energy values respecfully for plotting
particle electron(pos,2.1,0.5,it);	// Don't care about the arguments, just to initialize
electron.wavefunction(); 		// Initializes the variables in electron
double s;

//electron.generate_metropolis(metro_step);
observable energy_part(electron);

for (int j=0 ; j <= N+1; j++)
{	
	s=s_min+s_step*j;
	out = minimize(energy_part, beta_0, beta_1, eps,s); 
	dat1[j]=s;
	dat2[j]=out[1];
	fprintf(stdout,"%2.1f percent\n", 100*(float)j/(float)(N+1));
}


//////////////////    RESULTS     //////////////////////

mydata testclass(dat1, dat2, N);
testclass.plot2d();
return 0;
}

