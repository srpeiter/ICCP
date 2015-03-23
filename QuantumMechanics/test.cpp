#include "allheaders.h"

using namespace std;

int main(void)

{
////////////////  PARAMETERS   /////////////////

double metro_step = 0.1;
double s_min = 0.1;
double s_max = 3;
double s_step = 0.01;
int it = 100000; // Number of iterations in the Monte Carlo 
double beta_0 = 8.8; // Begining values for the minimiztion algorithm
double beta_1 = 9;
double conv_indicator = 0.1; // The minimization algorithm stops when the difference between two successive values is lower
double eps = 0.1*beta_0; // Speed of convergeance in minimization

////////////////   ALGORITHM   /////////////////

int N = (int)((s_max - s_min)/s_step); // Number of values of s
double pos[2][3]= {{1,1,1},{1,1,1}}; // Initial position of electrons
double *out; // Used to store the minimization results: the minimal energy and the beta at which this is achieved
double dat1[N], dat2[N]; // Stores the s values and the energy values respecfully for plotting
particle electron(pos,2.1,0.5,it);	// Don't care about the arguments, just to initialize
electron.wavefunction(); // Initializes the variables in electron

electron.generate_metropolis(metro_step);
observable energy_part(electron);

for (int j=0 ; j <= N+1; j++)
{
	out = minimize(energy_part, beta_0, beta_1, eps,(s_min+s_step*j)); 
	dat1[j]=(s_min+s_step*j);
	dat2[j]=out[1];
}


//////////////////    RESULTS     //////////////////////

mydata testclass(dat1, dat2, N);
testclass.plot2d();
return 0;
}

