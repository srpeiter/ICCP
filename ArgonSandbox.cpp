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

int N;			// Number of particles
int Run; 		// Number of simulations rounds
int run;		// Round number
int runtemp = 0;	// Used for the display function 
double dt;		// discrete time step
const int dim = 3;	// Number of dimension
double EK;		// Total kinetic energy
double EP;		// Total potential energy
double Boxlength =10; 	// The length of the box
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
void Initiate_Velocity();
void Update_position();
void Update_velocity();
void Update_force();
void Update_boundaries();
void Display();
void Calculate_Potential_Energy();
void Calculate_Kinetic_Energy();
void Energy_Correction();


/*
_________________________________________________________

 -------------- Random number generator ---------------- 
_________________________________________________________

*/

default_random_engine generator;
uniform_real_distribution<double> unif_dist(0.0,1.0);

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
if (argc < 3){
printf("Error: Give the number of particles N followed by the number of runs. \nExample: './ArgonSandbox 9 10 ', meaning 9 particles and 10 runs 1\n");}
else{
N = atoi(argv[1]);
Run = atoi(argv[2]); 
dt = 0.0001;
}
}

// ________ Initiate particles position ________ \\

void Initiate_Position_Simple(){
for (int n=0; n<N; ++n){
pos[0][n] = 1.2 * n;	// Let's start easy
pos[1][n] = 1;
pos[2][n] = 1;
}
}

void Initiate_Position_Cubic(){
int N3 = round(pow(N,0.3333));
int n = 0;
for (int x = 0; x < N3; ++x){
for (int y = 0; y < N3; ++y){
for (int z = 0; z < N3; ++z){
pos[0][n] = 1*x;	
pos[1][n] = 1*y;
pos[2][n] = 1*z;
n += 1;
} } }
}

// ________ Initiate particles velocity ________ \\

void Initiate_Velocity(){
for (int n=0; n<N; ++n){
vel[0][n] = 0;
vel[1][n] = 0;
vel[2][n] = 0;
}
}

/*
_________________________________________________________

 --------------- Updating functions -------------------- 
_________________________________________________________

*/

// ________ Update particles position ________ \\

void Update_position(){
for (int n = 0; n < N; ++n){
for (int i = 0; i < dim; ++i){		
pos[i][n] = pos[i][n] + vel[i][n] *dt;	// x(i) = x(i-1) + v(i)*dt
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


void Update_force(){
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
int Shift_X = -1;
int Shift_Y = -1;
int Shift_Z = -1;
for (int b = 0; b << 27; ++b){ // Box 1 till 27
if (b = 3)
cout << " test" << endl;
for (int n = 0; n << N; ++n){
POS[0][n+(b * N)] = pos[0][n] + Boxlength * Shift_X; 
POS[1][n+(b * N)] = pos[1][n] + Boxlength * Shift_Y; 
POS[2][n+(b * N)] = pos[2][n] + Boxlength * Shift_Z; 
} }
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
for (int n = 0; n < N; ++n){
for (int d = 0; d < dim; ++d)
EK += 0.5 * pow(vel[d][n],2);	 	// Calculate kinet energy /particle/dim
}
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
} }
if (EnergyDisp > 0){
cout << "EK: " << EK << " EP: " << EP << " EP + EK: " << EP+EK << endl;
}
}
runtemp += 1;
}

/*
_________________________________________________________

 -------------------- Main file ------------------------ 
_________________________________________________________

*/

int main(int argc, char*argv[]){
Parameters(argc, argv);
Make_array(argc, argv);
Initiate_Position_Cubic();
Initiate_Velocity();
for (run = 0; run<Run; ++run){
Display();
Update_force();
Update_velocity();
Update_position();
Calculate_Kinetic_Energy();
Calculate_Potential_Energy();
Energy_Correction();
}

}
