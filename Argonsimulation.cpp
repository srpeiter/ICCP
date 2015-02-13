#include<iostream>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<random>
// #include"dislin.h"

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
int Run = 10;		// Number of runs
float dt = 0.01;		// Discrete time unit
const int dim =3; // read-only value
double *pos[dim];	// position array x,y,z
double *old_pos[dim];
double *vel[dim];  // velocity array vx, vy, vz;
double *Force[dim];



const int cell_number=5;  // test: cellsize should be the cuttoff r_c of 
// LJ-potential, It's important since it eases our computation
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
void setup_cell_link();
void force_calculate();
void brute_force_calculate();
void Update_velocity();
void Update_position();
void Calculate_Energy();
void Calculate_Individual_Kinetic_Energies();
void Calculate_Individual_Potential_Energies();

void Make_array(int argc, char* argv[] )
{
int i;
if (argc!=2) // two inputs , which are executable and number of particles N
	printf("Error: Give the number of particles N");

else {
	boxlength=atoi(argv[1]); // extract number from command line, box                                      // starting from 0
        N=pow(boxlength,3);
	link_list =new int[N];

	//allocate memory: not fully dynamic, because rows are on stack thus 
	//fixed and columns are on heap, which is dynamic; So partially
	//dynamic
	for (i = 0; i< dim; i++){
	pos[i]= new double[N] (); 
	vel[i]= new double[N] (); 
	Force[i]= new double[N] ();
	old_pos[i] = new double[N] ();
	} 
}	
 
if (pos == NULL | vel==NULL | link_list==NULL | Force==NULL)
printf("Error:Memory couldn't be allocated");


//memset(pos,0,N*3*sizeof(pos[0][0]));
//memset(vel,0,N*3*sizeof(vel[0][0]));
//memset(Force,0,N*3*sizeof(Force[0][0]));
memset(link_list,0,sizeof(link_list));


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
int m=0,n=0;
double x,y,lx,ly; //fcc making variable

int stepsize=boxlength;
// making random number generator
// Using a normal distribution



	for (int i = 0; i < stepsize; i++)
		for(int j=0;j < stepsize; j++)
			for( int k=0; k < stepsize;k++)

			{  if (n>N)
				return;

		// making fcc lattice	
			
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

		//cout << "v_mean"<< " "<< v2_sum<< endl;
		// random number generator for initial velocity determination



}





// Step 2. We compute the forces between particles with the link list method
// this force computing step is most CPU-time consuming

//function to make cell list and link


void setup_cell_link()
{

int cellx,celly, cellz;  // cell index in HEAD

//initialize cell_matrix, thus positions of particles in 3D box

memset(HEAD,0,sizeof(HEAD[0][0][0]*pow(cell_number,3)));

 for (int i = 0; i < N; ++i) 
 {
 cellx= floor(cell_number*pos[0][i]/(boxlength*gridsize));
 celly= floor(cell_number*pos[1][i]/(boxlength*gridsize));
 cellz= floor(cell_number*pos[2][i]/(boxlength*gridsize));

link_list[i]=HEAD[cellx][celly][cellz];
//cout << link_list[i]<< endl;
HEAD[cellx][celly][cellz]=i;

}
}



// When Cell_link mastrices are set up, we now calculate the force
// on each particles



void force_calculate()
{

/* we start by locating the paricles in the box
 * and then look for other particles in that box and
 * particles in the neighbouring box (27 neighbours to be correct)
 * and calcute the force based on the lennard and jones Potential
 */

// first locate particles in box.
int cellx,celly,cellz;
// array to store if the is the period whether the periodicity of hte Boundary
double postemp[3];
double dist[3], radius;

int ii, jj,kk, id;

for (int m=0 ; m < N; ++m)
{// cout << "m "<< m<< endl;
 cellx= floor(cell_number*pos[0][m]/((boxlength)*gridsize));
 celly= floor(cell_number*pos[1][m]/((boxlength)*gridsize));
 cellz= floor(cell_number*pos[2][m]/((boxlength)*gridsize));
  

  for(int i=cellx-1; i <= cellx+1 ; i++)
	for(int j=celly-1; j <=celly+1 ;j++)
		for(int k=cellz-1; k <= cellz+1; k++)
		{
		int periodic[3]={0,0,0};
               
		// identify edges left/right and up/Down
		// periodic : 0 is no boundary, 1 is lower boundary and 2 
		// is upper Boundary
		ii=i;// cout <<"ii "<<ii<<endl;
		jj=j;// cout <<"jj "<<jj<<endl;

		kk=k; //cout <<"kk "<<kk<<endl;

		if(ii < 0)
		{
		periodic[0]=1;
		ii+= cell_number;
		}
		else if(ii>cell_number-1)
		{
		periodic[0]=2;
		ii-=cell_number;
		}
		
		if(jj < 0)
		{
		periodic[1]=1;
		jj+= cell_number;
		}
		else if(jj>cell_number-1)
		{
		periodic[1]=2;
		jj-=cell_number;
		}
		
		if(kk < 0)
		{
		periodic[2]=1;
		kk+= cell_number;
		}
		else if(kk>cell_number-1)
		{
		periodic[2]=2;
		kk-=cell_number;
		}
	
// now locate particle in box;

id = HEAD[ii][jj][kk];

// search for all particles in box;
do {

if (id == m)
id = link_list[id]; // skip if it is the same particles

// calculate distance based on where they are located in box ex. boundaries,
// middle
for (int d=0; d < dim; ++d)
{
if (periodic[d]==0)
	postemp[d]=pos[d][id];
else if (periodic[d]==1)
	postemp[d]=pos[d][id] +(boxlength-1) ;
else if (periodic[d]==2)
	postemp[d]=pos[d][id] - (boxlength-1);
dist[d]=pos[d][m]-postemp[d];
}
// calculate distance
radius= pow(dist[0],2) +  pow(dist[1],2) + pow(dist[2],2) ;

// now calculate force through lennard jones Potential
for (int d=0; d < dim; ++d)
Force[d][m] += (-24* (dist[d])*(2*pow(radius,-7)+ pow(radius,-4)));



id= link_list[id];
//cout << "id "<< id << endl;
} while (id !=0);// go till the end of the link_list	

}	
 
}
}

void brute_force_calculate(){
/* This file uses brute force (without optimalization) to calculate the  force between particles. This will be used to test the more optimal force calculation for their accuracy.
*/
int n,m;
double dist[3],radius;

for (n = 0; n < N; ++n){
for (int d=0; d < dim; ++d)
Force[d][n] = 0;	// Reset the force

for (m=0; m < N; ++m){
// Skip if the loop is dealing with the same particle
if (n != m){

// Distance between particle n and all other partles m
for (int d=0; d < dim; ++d)
{
dist[d] = pos[d][m]-pos[d][n]; 
}

// Calculate the distance (squared) between two particles
radius = pow(pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2),3);

// Calculate the force according to the lennard jones Potential
for (int d=0; d < dim; ++d)
{
Force[d][n] += 6 * (2*pow(radius,-2) - pow(radius,-1))*dist[d];
} 

} } } }

// Step 3. Update position and velocity 
//

// ------ Update the velocity ------ \\

void Update_velocity(){
int i,n;
for (n = 0; n < N; ++n){
for (i = 0; i < dim; ++i) {
vel[i][n] = vel[i][n] - Force[i][n] * dt; // Every time dt, the particle gets an 						  // increase in velocity[x,y,z]:force * dt
} } 
}

// ------ Update position ------ \\

void Update_position(){
int i,n;
for (n = 0; n < N; ++n){
for (i = 0; i < dim; ++i) {
old_pos[i][n] = pos[i][n]; 			// save previous position in new matrix
pos[i][n] = pos[i][n] + vel[i][n] *dt;
}}
}

// Step 4. Characterize system
// The system has certain quantities that we would like to extract. For example the Kinetic and Potential energy of the system.

// ------ Total Energy ------ \\

void Calculate_Energy(){
double Ek = 0;
double Ep = 0; 
double dist[3];
double radius;

//  Kinetical energy (v^2/2) \\

for (int n = 0; n < N; ++n)
Ek += pow(vel[0][n],2)+pow(vel[1][n],2)+pow(vel[2][n],2);
Ek = 0.5 * Ek;

//  Potential energy (Lennard jones potential) \\

for (int n = 0; n < N; ++n){
for (int m=0; m < N; ++m){
if (n != m){ // Skip if the loop is dealing with the same particle
// Distance between particle n and all other partles m
for (int d=0; d < dim; ++d)
dist[d] = pos[d][m]-pos[d][n]; 
// Calculate the distance (squared) between two particles
radius = pow(pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2),3);
// Calculate the potential energy according to the lennard jones Potential
Ep += ((pow(radius,-2) - pow(radius,-1)));
}}}
cout << Ep + Ek << endl;
}

// ------ Kinetic energy of every particle ------ \\

void Calculate_Individual_Kinetic_Energies(){
double Ek[N],EK=0;
for (int n=0; n<N; ++n)
Ek[n] = 0.5 * pow(vel[0][n],2)+pow(vel[1][n],2)+pow(vel[2][n],2);
for (int n=0; n < N ;n++)
EK += Ek[n];
cout << "EK " << EK << endl;
}

// ------ Potential energy of every particle ------ \\

void Calculate_Individual_Potential_Energies(){
double radius;
double Ep[N],EP=0;
double dist[dim];
for (int n=0; n<N; ++n){
Ep[n]=0;
for (int m=0; m<N; ++m){
if (n != m){
radius = 0;
for (int d=0; d<dim; ++d){
dist[d] = pow(pos[d][n]-pos[d][m],2); 
radius += dist[d];
}
Ep[n] += pow(radius,-6) - pow(radius,-3);
} } 
EP +=Ep[n];
}
cout << "EP " << EP <<endl;
}



int main(int argc, char* argv[])
{
int i;


// Initiate simulation
 Make_array(argc,argv);
 Initialization();
 //setup_cell_link();

/* Run simulation
The order is important!
The force is dependend on the position, the velocity is dependend on the force and the position is dependend on the velocity.
So: First calculate the force according to the old positions. followed by the update of the position and velocity after that. 
Note: I'm not sure if the previous force or the new force should be used to update teh velocity.
*/

for (i = 0; i < Run ; i++){
//force_calculate();

Update_position();	
Update_velocity(); 
brute_force_calculate();
//Calculate_Energy();
Calculate_Individual_Kinetic_Energies();
Calculate_Individual_Potential_Energies();

//int sim=0;
//sim++;
//if (sim>100){
//sim = sim -100;
//cout << i << " " << Force[0][2] << " " << Force[1][2] << " " << Force[2][2] << endl;
//cout << i << " " << pos[0][2] << " " << pos[1][2] << " " << pos[2][2] << endl;
//}



}





// Debugging tools:
//for (int n=0; n < N ;n++)
//cout << pos[0][n] << " "<< pos[1][n]<< " "<< pos[2][n]<<endl;
//cout << Force[0][n] << " "<< Force[1][n]<< " "<< Force[2][n]<<endl;


for (int i=0; i<3; i++)
{delete [] pos[i];
delete [] vel[i];
delete [] Force[i];
}
delete [] link_list;
}
