#include<iostream>
#include<math.h>
#include<stdio.h>
#include<random>

using namespace std;

/*In this code we simualte motion of gas particles that interact 
 * through the Lennard-Jones Potential.
 * We implement this code in several steps:
 * 1. Initialization : position and velocity
 * 2. Motion: This computed based on Newton equations
 * 3. Boundary conditions : Periodic and/or fixed
 * 4. Obtaining phase diagrams
 */

// Global position variable
int N;  		// number of particles 
const int dim =3; // read-only value
double *pos[dim];	// position array x,y,z
double *vel[dim];  // velocity array vx, vy, vz;

// function to determine the number of particles during runtime
// and then to make the array pos and vel (allocate/reserve memory);

void Make_array();
void Initialization();



void Make_array(int argc, char* argv[] )
{
int i;
if (argc!=2) // two inputs , which are executable and number of particles N
	printf("Error: Give the number of particles N");

else {
	N=atoi(argv[1]); // extract number from command line

	//allocate memory: not fully dynamic, because rows are on stack thus 
	//fixed and columns are on heap, which is dynamic; So partially
	//dynamic
	for (i = 0; i< dim; i++){
	pos[i]= new double[N];
	vel[i]= new double[N];
	
	} 	
if (pos == NULL | vel==NULL)
printf("Error:Memory couldn't be allocated");
}
}


//Step 1: Initialization
// This function takes Global variable: position and velocity array
// and sort this in a cubic box (fcc structure) and gives them a velocity based
// on the maxwell velocity distribution.
// You only need to define the length of the array, which is equal to the 
// number of particles N.
void Initialization()
{
 // Particles initial positions are chosen such, that there aren't any overlap
 // So place the particles on a cubic lattice.
 // box has length of 10 and for example if N=1000, we put each particles
 // 1 unit length from each other
 // this routine is done once, so we afford to use some for loops

int i,j,k,n=0, boxlength=10;
int stepsize = ceil(pow(N,1/3));
double gridsize= boxlength/stepsize;	// gridsize position where the particles
					// should be placed;
// making random number generator
// Using a normal distribution

default_random_engine generator;
normal_distribution<double> maxwell(0,1);

	for (i = 1; i < stepsize; i++)
		for(j=1;j < stepsize; j++)
			for(k=1; k < stepsize;k++)

			{  if (n>N)
				return;
			//position initialization
			pos[0][n]=(3 * i + j + k) * gridsize;
			pos[1][n]=(3 * j + h) * gridsize;
			pos[2][n]=(3 * h) * gridsize;

			// velocity initialization
			vel[0][n]=maxwell(generator);
			vel[1][n]=maxwell(generator);
			vel[2][n]=maxwell(generator);
			
			++n;
			
			
		
	}
}


// Step 2. We compute the forces between particles with the link list method
//
int main(int argc, char* argv[])
{

 Make_array(argc,argv);
 Initialization();


cout<<"Whats up"<< endl;

}
