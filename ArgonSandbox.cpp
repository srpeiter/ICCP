#include<iostream>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<random>

using namespace std;

/* 
This is a kind or remake of Argonsimulation.cpp.
Out of hopelessness, I created a new file which is more simple such that I can investigate the problem more isolated.
*/

/* 
NOTES:
The function Distance is not pralocated correctly
The functions Distance and force do not work yet
*/


/*
_________________________________________________________

 ------------ Global variable initiation --------------- 
_________________________________________________________

*/

int N;			// Number of particles
int Run; 		// Number of simulations rounds
double dt;		// discrete time step
const int dim = 3;	// Number of dimensions
double *pos[dim];	// Position array x,y,z
double *vel[dim];	// Velocity array x,y,z
double *forc[dim];	// Force array x,y,z

double *dist[dim];	// A temperary matrix used for distance calculations
//double Dist2[N];	// A temperary vector also used for distane calculations

/*
_________________________________________________________

 -------------- Function initialization ---------------- 
_________________________________________________________

*/

void Make_array();
void Parameters();
void Initiate_Position();
void Initiate_Velocity();
void Display();
void Distance();

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
if (argc > 2){
N=atoi(argv[1]); 
for (i = 0; i< dim; i++){
pos[i]= new double[N] (); 
vel[i]= new double[N] (); 
forc[i]= new double[N] ();
dist[i] = new double[N] ();
} } 
}

// ________ Process parameters ________ \\

void Parameters(int argc, char* argv[]){
if (argc < 2){
printf("Error: Give the number of particles N followed by the number of runs. \nExample: './ArgonSandbox 9 10', meaning 9 particles and 10 runs\n");}
else 
Run = atoi(argv[2]); 
}

// ________ Initiate particles position ________ \\

void Initiate_Position(){
for (int n=0; n<N; ++n){
pos[0][n] = n;	// Let's start easy
pos[1][n] = 0;
pos[2][n] = 0;
}
}

// ________ Initiate particles velocity ________ \\

void Initiate_Velocity(){
for (int n=0; n<N; ++n){
vel[0][n] = 0;
vel[1][n] = n;
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
pos[i][n] = pos[i][n] + vel[i][n] *dt;
}}
}

// ________ Update particles velocity ________ \\


void Update_velocity(){
for (int n = 0; n < N; ++n){
for (int i = 0; i < dim; ++i){		
vel[i][n] = vel[i][n] + forc[i][n] *dt;
}}
}

// ________ Distance ________ \\


void Distance(int& n){
cout << "test" << endl;
/*
for (int m = 0; m < N; ++m){
for (int d = 0; d < dim; ++d) {		
dist[d][m] = pos[d][m] - pos[d][n];
}
Dist2[m] = pow(dist[0][m],2) + pow(dist[1][m],2) + pow(dist[2][m],2)
}
*/
}

// ________ Update particles force ________ \\


void Update_force(){
int i,n;
for (n = 0; n < N; ++n){
for (i = 0; i < dim; ++i) {
Distance(n);
		
//forc[i][n] = forc[i][n] + forc[i][n] *dt;
}}
}



/*
_________________________________________________________

 --------------- Extra functions -------------------- 
_________________________________________________________

*/

// ________ Display values on screen ________ \\

void Display(){
int n;
for (n = 0; n < N; ++n)		
cout << pos[1][n] << " " << endl ;
}


/*
_________________________________________________________

 -------------------- Main file ------------------------ 
_________________________________________________________

*/

int main(int argc, char*argv[]){
Parameters(argc, argv);
cout << "Test" ;
Make_array(argc, argv);
Initiate_Position();
Initiate_Velocity();
for (int run = 0; run<Run; ++run){
dt = 0.01;
Update_position();
Display();
}

}
