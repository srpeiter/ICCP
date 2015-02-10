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
double gridsize= 1;   		// gridsize position where the particles
				// should be placed;
const int dim =3; // read-only value
double *pos[dim];	// position array x,y,z
double *vel[dim];  // velocity array vx, vy, vz;



const int cell_number=10;  // test: cellsize should be the cuttoff r_c of 
		// LJ-potential, 
// Cell and link matrices associated with force calculation: short range 
// interactions

int HEAD[cell_number][cell_number][cell_number];
int *link_list;

// random number generator for initial velocity determination
default_random_engine generator;
uniform_real_distribution<double> unif_dist(0.0,1.0);

// function to determine the number of particles during runtime
// and then to make the array pos and vel (allocate/reserve memory);
//
void Make_array();
double box_muller_gauss();
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

//polar form of box-muller for generating Maxwell distribution
// mean 0 and std 1;
double  box_muller_gauss()
{
double ran_number, x1, x2, w;


do {

x1= 2* unif_dist(generator) -1;
x2= 2* unif_dist(generator) -1;
w = x1*x1 + x2*x2;
} while (w >=1);

return ran_number= x1* (sqrt((-2 * log(w))/w));

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
double Temp=5, sf; // sf is scaling vactor;
double v_sum[3] = {0,0,0};
double v2_sum=0.0;
int m=0,n=0,x,y,lx,ly;

int stepsize=boxlength;
// making random number generator
// Using a normal distribution



	for (int i = 0; i < stepsize; i++)
		for(int j=0;j < stepsize; j++)
			for( int k=0; k < stepsize;k++)

			{  if (n>N)
				return;

			
			
			 //position initialization
			x = (3*i + j + k)/3;	//x position 
			y = (3*j + k)/3;	//y position
						//z = k
			   //subtract step size if x or y > step size
			lx = stepsize * floor(x/stepsize); 	
			ly = stepsize * floor(y/stepsize);
			pos[0][n]=(x-lx)*gridsize;
			pos[1][n]=(y-ly)*gridsize;
			pos[2][n]=(k)*gridsize;
			++n;
			}

	// velocity Initialization
			for ( m=0; m < N ;++m){
			vel[0][m]=box_muller_gauss();
			vel[1][m]=box_muller_gauss();
			vel[2][m]=box_muller_gauss();
			v_sum[0]+= vel[0][m];
			v_sum[1]+= vel[1][m];
			v_sum[2]+= vel[2][m];
	v2_sum += pow(vel[0][m],2)+ pow(vel[1][m],2)+ pow(vel[2][m],2);
			}

			v_sum[0]=v_sum[0]/N;
			v_sum[1]=v_sum[1]/N;
			v_sum[2]=v_sum[2]/N;
			v2_sum= v2_sum/N;       // rescaling velocity
			sf= sqrt(3*Temp/v2_sum) ;// rescaling factor

		m=0;
	        v2_sum=0.0;

		for (m = 0; m < N; m++) {
		vel[0][m]=(vel[0][m]-v_sum[0])*sf;
		vel[1][m]=(vel[1][m]-v_sum[1])*sf;
		vel[2][m]=(vel[2][m]-v_sum[2])*sf;
	v2_sum += pow(vel[0][m],2)+ pow(vel[1][m],2)+ pow(vel[2][m],2);

		} 

		cout << "v_mean"<< " "<< v2_sum<< endl;
		// random number generator for initial velocity determination



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
 cellx= floor(cell_number*pos[0][i]/(boxlength*gridsize));
 celly= floor(cell_number*pos[1][i]/(boxlength*gridsize));
 cellz= floor(cell_number*pos[2][i]/(boxlength*gridsize));

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
cout << vel[0][n] << " "<< vel[1][n]<< " "<< vel[2][n]<<endl;

 cell_link();


for (int i=0; i<3; i++)
{delete [] pos[i];
delete [] vel[i];
}
delete [] link_list;
}
