#include<iostream>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<random>

using namespace std;

/* 
This is a kind or remake of Argonsimulation.cpp.
Out of hopelessness, I created a new file that startes simpler such that I can investigate the problem more isolated.
*/

/* 
NOTES:
Life is awesome
*/


/*
_________________________________________________________

 ------------ Global variable initiation --------------- 
_________________________________________________________

*/
double a;	   	// Length of unit FCC unit cell
int N;			// Number of particles
int cub_num;		// 
int Run; 		// Number of simulations rounds
int run;		// Round number
int runtemp = 0;	// Used for the display function 
double dt;		// discrete time step
double density;		// The denisty of the box
const int dim = 3;	// Number of dimension
double EK;		// Total kinetic energy
double EP;		// Total potential energy
double boxlength; 	// The length of the box
double *pos[dim];	// Position array X,Y,Z
double *POS[dim];	// Position array X,Y,X with boundaries included
double *vel[dim];	// Velocity array Vx,Vy,Vz
double *force[dim];	// Force array Fx,Fy,Fz

double *dist[dim];	// A temperary matrix used for distance calculations


/*
_________________________________________________________

 -------------- Function initialization ---------------- 
_________________________________________________________

*/

void Make_array();
void Parameters();
void Initiate_Position_Simple();
void Initiate_Position_Cubic();
void Initiate_Position_FCC();
void Initiate_Velocity();
void Update_position();
void Update_velocity();
void Update_position_with_boundaries();
void Update_force_brute_force();
void Update_force_with_boundaries_brute_force();
void Update_boundaries();
void Display();
void Matrix_row_checker();
void Calculate_Potential_Energy();
void Calculate_Kinetic_Energy();
void Energy_Correction();


/*
_________________________________________________________

 -------------- Random number generator ---------------- 
_________________________________________________________

*/

default_random_engine generator;
normal_distribution<double> maxwell(0.0,1.0);

/*
_________________________________________________________

 --------------- Initiating Functions ------------------ 
_________________________________________________________

*/

// ________ Initiate matrices ________ \\

void Make_array(int argc, char* argv[])
{
int i;
// Matrix preallocating
for (i = 0; i< dim; i++){
pos[i]= new double[N] (); 
vel[i]= new double[N] (); 
force[i]= new double[N] ();
dist[i] = new double[N] ();
POS[i] = new double[27 * N] (); // The positions will be dublicated all around the box
} 
// Vector preallocating
}

// ________ Process parameters ________ \\

void Parameters(int argc, char* argv[]){
if (argc !=5){
printf("Error: Give the number of particles N followed by the denisty, number of runs and time steps. \nExample: './ArgonSandbox 9 1 10 0.0001', meaning 9 particles and 10 runs with density 1. Taking time steps dt.\n");}
else{

// Inputes
N = atoi(argv[1]);
density = atof(argv[2]);
Run = atoi(argv[3]); 
dt = atof(argv[4]);

// Others and calculated values
boxlength = pow((N/density),0.3333);
cub_num = ceil(pow(N,0.333));
a=boxlength/cub_num;
}
}

// ________ Initiate particles position ________ \\

void Initiate_Position_Simple(){
for (int n=0; n<N; ++n){
pos[0][n] = 1.2 * n;	
pos[1][n] = 1;
pos[2][n] = 1;
}
}

void Initiate_Position_Cubic(){
int n = 0;
for (int x = 0; x < cub_num; ++x)
for (int y = 0; y < cub_num; ++y)
for (int z = 0; z < cub_num; ++z)
if (n<N){ 			// Stop when the number of particles is reached
pos[0][n] = x * a;	
pos[1][n] = y * a;
pos[2][n] = z * a;
n++;
}
}

void Initiate_Position_FCC(){
double unitcell[3][4]={{0, 0.5*a, 0.5*a,0},{0,0,0.5*a ,0.5*a},{0,0.5*a,0,0.5*a}};
int n = 0;
for (int x = 0; x < cub_num; ++x)
for (int y = 0; y < cub_num; ++y)
for (int z = 0; z < cub_num; ++z)
for (int j = 0; j < 4; ++j)
if (n < N){ 			// Stop when the number of particles is reached
pos[0][n] = unitcell[0][j] + x * a;	
pos[1][n] = unitcell[1][j] + y * a;
pos[2][n] = unitcell[2][j] + z * a;
n++;
}
}


// ________ Initiate particles velocity ________ \\

void Initiate_Velocity(){
for (int n=0; n<N; ++n){
vel[0][n] = maxwell(generator);
vel[1][n] = maxwell(generator);
vel[2][n] = maxwell(generator);
}
}

/*
_________________________________________________________

 --------------- Updating functions -------------------- 
_________________________________________________________

*/

// ________ Update particles position ________ \\

void Update_position(){

// basic function

for (int n = 0; n < N; ++n){
for (int i = 0; i < dim; ++i){		
pos[i][n] += pos[i][n] + vel[i][n] *dt;	// x(i) = x(i-1) + v(i)*dt
}}
}

void Update_position_with_boundaries(){

// This function puts particles back in the box when they would otherwise leave the box

for (int n = 0; n < N; ++n){
for (int i = 0; i < dim; ++i){		
pos[i][n] += vel[i][n] *dt;	// x(i) = x(i-1) + v(i)*dt
pos[i][n] += -floor(pos[i][n]/boxlength); // subtracts (-1, 0 or 1) * boxlength
}}
}

// ________ Update particles velocity ________ \\


void Update_velocity(){
for (int n = 0; n < N; ++n){
for (int i = 0; i < dim; ++i){		
vel[i][n] = vel[i][n] + force[i][n] *dt;// v(i) = v(i-1) + f(i)*dt
}}
}

// ________ Update particles force ________ \\


void Update_force_with_boundaries_brute_force(){

// This function calculates the force between every particle for each particle including the particles due to the boundary conditions (except itself) 

double Dist2;
int d;
for (int n = 0; n < N; ++n){
for (d = 0; d < dim; ++d)
force[d][n] = 0; 			// Reset the force
for (int m = 0; m < 27 * N; ++m){	// Including boundaries, 27 x N particles
//if (13 * N + n != m){			// Only calculate between different 
					// particles (mind that we are dealing with 
					// comparing with the boundary matrix 
					// (middle box needed)
Dist2 = 0;				// Reset "Distance squared"
for (d = 0; d < dim; ++d){
dist[d][m] = POS[d][m] - pos[d][n];	// Calculate the distance in x,y or z
Dist2 += pow(dist[d][m],2);		// Calculate the distance squared
}
if (Dist2 != 0){
for (d = 0; d < dim; ++d)		// Calculate the force
force[d][n] += 24 * (-2*pow(Dist2,-7) + pow(Dist2,-4)) * dist[d][m];
} } }
}

void Update_force_brute_force(){

// This function calculates the force between every particle for each particle (except itself) 

double Dist2;
int d;
for (int n = 0; n < N; ++n){
for (d = 0; d < dim; ++d)
force[d][n] = 0; 			// Reset the force
for (int m = 0; m < N; ++m){
if (n !=m){				// Only calculate between different 
					// particles
Dist2 = 0;				// Reset "Distance squared"
for (d = 0; d < dim; ++d){
dist[d][m] = pos[d][m] - pos[d][n];	// Calculate the distance in x,y or z
Dist2 += pow(dist[d][m],2);		// Calculate the distance squared
}
for (d = 0; d < dim; ++d)		// Calculate the force
force[d][n] += 24 * (-2*pow(Dist2,-7) + pow(Dist2,-4)) * dist[d][m];
} } }
}

// ________ Update boundaries ________ \\

// copy the particles in the box all around the middle box
// Note: under construction

void Update_boundaries(){
int b = 0;
for (int x = -1; x < 2; ++x){
for (int y = -1; y < 2; ++y){
for (int z = -1; z < 2; ++z){
for (int n = 0; n < N; ++n){
POS[0][n+(b * N)] = pos[0][n] + boxlength * x; 
POS[1][n+(b * N)] = pos[1][n] + boxlength * y; 
POS[2][n+(b * N)] = pos[2][n] + boxlength * z; 
} 
b++;
} } } 
}


/*
_________________________________________________________

 ------------ Energy related functions ----------------- 
_________________________________________________________

*/

// ________ Calculate Potential Energy ________ \\

// Note: It is more optimal to calculate the force and potential energy in one function

void Calculate_Potential_Energy(){
double Dist2;
EP = 0; 				// Reset total potential energy parameter
for (int n = 0; n < N; ++n){
for (int m = n+1; m < N; ++m){		// Do not count double the potential energy
Dist2 =0;				// Reset distance squared
for (int d = 0; d < dim; ++d){
dist[d][m] = pos[d][m] - pos[d][n];	// Distance between particles in one 
					// direction
Dist2 += pow(dist[d][m],2); 		// Distance between particles
}
EP += 4 * (pow(Dist2,-6)-pow(Dist2,-3)); // Total potential energy
} }
}

// ________ Calculate Kinetic Energy ________ \\

void Calculate_Kinetic_Energy(){
EK = 0;      				// Reset total kinetic energy parameter
for (int n = 0; n < N; ++n)
for (int d = 0; d < dim; ++d)
EK += 0.5 * pow(vel[d][n],2);	 	// Calculate kinet energy /particle/dim
}

// ________ Energy Correction ________ \\

void Energy_Correction(){

}

/*
_________________________________________________________

 --------------- Extra functions -------------------- 
_________________________________________________________

*/

// ________ Display values on screen ________ \\

void Display(){
int EnergyDisp = 1;		// Display energies (1/0)
int PVFDisp = 1;		// Display position, velocity, force (1/0)
if (runtemp > 100){		// only display every x runs
runtemp += -100;  		// reset counter
if (PVFDisp > 0){
for (int n = 0; n < N; ++n){
// cout << n << "." << run << " F: " << force[0][n] << " V: " << vel[0][n] << " X: " << pos[0][n] << endl;
// cout << n << "." << run << " X: " << pos[0][n] << " " << "Y: " << pos[1][n] << " " << "Z: " << pos[2][n] << " " << endl;
} 
for (int n = 0; n < 27 * N; ++n){
// cout << n << "." << run << " X: " << POS[0][n] << " " << "Y: " << POS[1][n] << " " << "Z: " << POS[2][n] << " " << endl;
} }
if (EnergyDisp > 0){
cout << "EK: " << EK << " EP: " << EP << " EP + EK: " << EP+EK << endl;
} }
runtemp += 1;
}

// ________ Check for matrix for rows with the same values ________ \\

void Matrix_row_checker(){
for (int n = 0; n < 27 * N; ++n)
for (int m = n + 1; m < 27 * N; ++m)
if (POS[0][n] == POS[0][m])
if (POS[1][n] == POS[1][m])
if (POS[2][n] == POS[2][m])
cout << n << " " << m << endl;
}

/*
_________________________________________________________

 -------------------- Main file ------------------------ 
_________________________________________________________

*/

int main(int argc, char*argv[]){
Parameters(argc, argv);
Make_array(argc, argv);
Initiate_Position_FCC();
Initiate_Velocity();
for (run = 0; run<Run; ++run){
Update_boundaries();
Display();
Update_force_with_boundaries_brute_force();
Update_velocity();
Update_position_with_boundaries();
Calculate_Kinetic_Energy();
Calculate_Potential_Energy();
Energy_Correction();
}

}
