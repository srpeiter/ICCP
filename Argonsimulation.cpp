#include<iostream>
#include<math.h>
#include<stdio.h>
#include<string.h>
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
int boxlength;
const int dim =3; // read-only value
double *pos[dim];	// position array x,y,z
double *vel[dim];  // velocity array vx, vy, vz;



const int cell_number=5;  // test: cellsize should be the cuttoff r_c of 
		// LJ-potential, 
// Cell and link matrices associated with force calculation: short range 
// interactions

int HEAD[cell_number][cell_number][cell_number];
int *link_list;

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
	boxlength=atoi(argv[1]); // extract number from command line
        N=pow(boxlength,3);
	link_list =new int[N];

	//allocate memory: not fully dynamic, because rows are on stack thus 
	//fixed and columns are on heap, which is dynamic; So partially
	//dynamic
	for (i = 0; i< dim; i++){
	pos[i]= new double[N];
	vel[i]= new double[N];
	
	} 	
if (pos == NULL | vel==NULL | link_list==NULL)
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

int i,j,k,n=0;

int stepsize=boxlength;
double gridsize= 0.5;	// gridsize position where the particles
					// should be placed;
// making random number generator
// Using a normal distribution

default_random_engine generator;
normal_distribution<double> maxwell(0,1);

	for (i = 0; i < stepsize; i++)
		for(j=0;j < stepsize; j++)
			for(k=0; k < stepsize;k++)

			{  if (n>N)
				return;
			//position initialization
			pos[0][n]=(i)*gridsize;
			pos[1][n]=(j)*gridsize;
			pos[2][n]=(k)*gridsize;
			
				// velocity initialization
			vel[0][n]=maxwell(generator);
			vel[1][n]=maxwell(generator);
			vel[2][n]=maxwell(generator);
			
			++n;
	
	}
}


// Step 2. We compute the forces between particles with the link list method
// this force computing step is most CPU-time consuming

//function to make cell list and link


void cell_link()
{

int cellx,celly, cellz;  // cell index in HEAD

//initialize cell_matrix, thus positions of particles in 3D box
memset(link_list,0,sizeof(link_list));
memset(HEAD,0,sizeof(HEAD[0][0][0]*pow(cell_number,3)));

 for (int i = 0; i < N; i++) 
 {
 cellx= floor(cell_number*pos[0][i]/(boxlength*0.5));
 celly= floor(cell_number*pos[1][i]/(boxlength*0.5));
 cellz= floor(cell_number*pos[2][i]/(boxlength*0.5));

link_list[i]=HEAD[cellx][celly][cellz];
cout << link_list[i]<< endl;
HEAD[cellx][celly][cellz]=i;


 }




}

int main(int argc, char* argv[])
{

 Make_array(argc,argv);
 Initialization();
for (int n=0; n < N ;n++)
cout << pos[0][n] << " "<< pos[1][n]<< " "<< pos[2][n]<<endl;

 cell_link();


for (int i=0; i<3; i++)
{delete [] pos[i];
delete [] vel[i];
}
delete [] link_list;
}
