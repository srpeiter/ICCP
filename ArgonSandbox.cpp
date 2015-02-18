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
double *pos[dim];	// Position array x,y,z
double *vel[dim];	// Velocity array x,y,z
double *force[dim];	// Force array x,y,z

double *dist[dim];	// A temperary matrix used for distance calculations


/*
_________________________________________________________

 -------------- Function initialization ---------------- 
_________________________________________________________

*/

void Make_array();
void Parameters();
void Initiate_Position();
void Initiate_Velocity();
void Update_position();
void Update_velocity();
void Update_force();
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
} 
// Vector preallocating
}

// ________ Process parameters ________ \\

void Parameters(int argc, char* argv[]){
if (argc < 2){
printf("Error: Give the number of particles N followed by the number of runs. \nExample: './ArgonSandbox 9 10 0.001', meaning 9 particles and 10 runs 1\n");}
else{
N = atoi(argv[1]);
Run = atoi(argv[2]); 
dt = 0.001;
}
}

// ________ Initiate particles position ________ \\

void Initiate_Position(){
for (int n=0; n<N; ++n){
pos[0][n] = 1.1 * n;	// Let's start easy
pos[1][n] = 0;
pos[2][n] = 0;
}
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
pos[i][n] = pos[i][n] + vel[i][n] *dt;		// x(i) = x(i-1) + v(i)*dt
}}
}

// ________ Update particles velocity ________ \\


void Update_velocity(){
for (int n = 0; n < N; ++n){
for (int i = 0; i < dim; ++i){		
vel[i][n] = vel[i][n] + force[i][n] *dt;	// v(i) = v(i-1) + f(i)*dt
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
for (int d = 0; d < dim; ++d)
for (int m = n+1; m < N; ++m){		// Do not count double the potential energy
Dist2 =0;				// Reset distance squared
for (int d = 0; d < dim; ++d){
dist[d][m] = pos[d][m] - pos[d][n];	// Distance between particles in one 
					// direction
Dist2 += pow(dist[d][m],2); 		// Distance between particles
}
EP += 4 * (pow(Dist2,-6) - pow(Dist2,-3)); 	// Total potential energy
}				
}
}

// ________ Calculate Kinetic Energy of each particle ________ \\

void Calculate_Kinetic_Energy(){
double Ek[N];
EK = 0;      				// Reset total kinetic energy parameter
for (int n = 0; n < N; ++n){
Ek[n] = 0; 				// Reset kinetic energy parameters
for (int d = 0; d < dim; ++d)
Ek[n] += 0.5 * pow(vel[n][d],2); 	// Calculate kinet energy /particle/dim
EK += Ek[n];
}
}

// ________ Calculate Kinetic Energy of each particle________ \\

void Energy_Correction(){

}

/*
_________________________________________________________

 --------------- Extra functions -------------------- 
_________________________________________________________

*/

// ________ Display values on screen ________ \\

void Display(){
if (runtemp > 10){
runtemp += -10;	
for (int n = 0; n < N; ++n){
cout << n << "." << run << " F: " << force[0][n] << " V: " << vel[0][n] << " X: " << pos[0][n] << endl;
}
cout << "EK: " << EK << " EP: " << EP << " EP + EK: " << EP+EK << endl;
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
Initiate_Position();
Initiate_Velocity();
for (run = 0; run<Run; ++run){
Update_force();
Update_velocity();
Update_position();
Calculate_Kinetic_Energy();
Calculate_Potential_Energy();
Energy_Correction();
Display();
}

}
