#include<iostream>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<random>

using namespace std;

/* 
NOTES:
Life is awesome
*/


/*
_________________________________________________________

 ------------ Global variable initiation --------------- 
_________________________________________________________

*/

double EK;		// Total kinetic energy
double EP;		// Total potential energy
double boxlength; 	// The length of the box
double dt;		// discrete time step
double density;		// The denisty of the box
double a;	   	// Length of unit FCC unit cell
double kb = 1.38066*pow(10,-23);	// Boltzmann constant
double Temp;		// Temperature

int N;			// Number of particles
int cub_num;		// Number of unit FCC boxes
int Run; 		// Number of simulations rounds
int run;		// Round number
int runtemp = -1;	// Used for the display function 
const int dim = 3;	// Number of dimension

double *pos[dim];	// Position array X,Y,Z
double *posb[dim];	// Position array X,Y,X with boundaries included
			// The size of this matrix is 27 * N
double *vel[dim];	// Velocity array Vx,Vy,Vz
double *force[dim];	// Force array Fx,Fy,Fz



/*
_________________________________________________________

 -------------- Function initialization ---------------- 
_________________________________________________________

*/

// Initiation Functions
void Make_array();
void Parameters();
void Initiate_Position_Simple();
void Initiate_Position_Cubic(); 			// Broken
void Initiate_Position_FCC();
void Initiate_Velocity();

// Updating functions
void Update_position();
void Update_velocity();
void Update_position_with_boundaries();
void Update_force_brute_force();
void Update_force_with_boundaries_brute_force();
void Update_boundaries();

// Energy related functions
void Calculate_Potential_Energy();
void Calculate_Potential_Energy_with_boundaries();
void Calculate_Kinetic_Energy();
void Energy_Correction(); 				// Empty

// Gas quantities
void Calculate_Temperature();
void Adjust_Temperature();

// Extra functions
void Clean_up();
void Display();
void Matrix_row_checker();
/*
_________________________________________________________

 -------------- Random number generator ---------------- 
_________________________________________________________

*/

default_random_engine generator;
normal_distribution<double> maxwell(0.0,1.0);

/*
_________________________________________________________

 --------------- Initiation Functions ------------------ 
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
posb[i] = new double[27 * N] (); // The positions will be duplicated all around the box
} 
// Vector preallocating
}

// ________ Process parameters ________ \\

void Parameters(int argc, char* argv[]){
if (argc !=5){
printf("Error: Give the number of particles N followed by the denisty, number of runs and time steps. \nExample: './ArgonSandbox 9 1 10 0.001', meaning 9 particles and 10 runs with density 1. Taking time steps 0.001.\n");}
else{

// Inputes
N = atoi(argv[1]);
density = atof(argv[2]);
Run = atoi(argv[3]); 
dt = atof(argv[4]);

// Others and calculated values
boxlength = pow((N/density),0.3333); 	// Lengths of the box
cub_num = ceil(pow((N/4),0.333));	// Number of boxes
a=boxlength/cub_num;			// Length of unity box
}
}

// ________ Initiate particles position ________ \\

// Sandbox style initialization

void Initiate_Position_Simple(){
for (int n=0; n<N; ++n){
pos[0][n] = 1.2 * n*a;	
pos[1][n] = 1;
pos[2][n] = 1;
}
}

// Cubic structure

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

// FCC structure

void Initiate_Position_FCC(){
double unitcell[3][4]={{0, 0.5*a, 0.5*a,0},{0,0,0.5*a ,0.5*a},{0,0.5*a,0,0.5*a}};
int n = 0;
for (int x = 0; x < cub_num; ++x)
for (int y = 0; y < cub_num; ++y)
for (int z = 0; z < cub_num; ++z)
for (int j = 0; j < 4; ++j)
if (n < N){ 			// Stop if the number of particles is reached
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


// basic function

void Update_position(){
for (int n = 0; n < N; ++n){
for (int i = 0; i < dim; ++i){		
pos[i][n] += vel[i][n] *dt;	// x(i) = x(i-1) + v(i)*dt
}}
}

// This function puts particles back in the box when they would 
// otherwise leave the box

void Update_position_with_boundaries(){
for (int n = 0; n < N; ++n){
for (int d = 0; d < dim; ++d){		
pos[d][n] += vel[d][n] *dt;	// x(i) = x(i-1) + v(i)*dt
pos[d][n] += -boxlength * floor(pos[d][n]/boxlength); // subtracts (-1, 0 or 1) * boxlength
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


// This function calculates the force between every particle for each 
// particle (except itself) 
// Note: "dist" as array is not necessary

void Update_force_brute_force(){
double Dist2;
double dist[dim];
int d;
for (int n = 0; n < N; ++n){
for (d = 0; d < dim; ++d)
force[d][n] = 0; 			// Reset the force
for (int m = 0; m < N; ++m){
if (n !=m){				// Only calculate between different 
					// particles
Dist2 = 0;				// Reset "Distance squared"
for (d = 0; d < dim; ++d){
dist[d] = pos[d][m] - pos[d][n];	// Calculate the distance in x,y or z
Dist2 += pow(dist[d],2);		// Calculate the distance squared
}
for (d = 0; d < dim; ++d){		// Calculate the force
force[d][n] += 24 * (-2*pow(Dist2,-7) + pow(Dist2,-4)) * dist[d];
}
} } }
}

// This function calculates the force between every particle for each 
// particle including the particles due to the boundary conditions 
// (except itself) 

void Update_force_with_boundaries_brute_force(){
double Dist2;
int d;
double dist[dim]; 
for (int n = 0; n < N; ++n){
for (d = 0; d < dim; ++d)
force[d][n] = 0; 			// Reset the force
for (int m = 0; m < 27 * N; ++m){	// Including boundaries, 27 x N particles
if (13 * N + n != m){			// Only calculate between different 
					// particles (mind that we are 
					// comparing with the boundary matrix 
					// (middle box needed)
Dist2 = 0;				// Reset "Distance squared"
for (d = 0; d < dim; ++d){
dist[d] = posb[d][m] - pos[d][n];	// Calculate the distance in x,y or z
Dist2 += pow(dist[d],2);		// Calculate the distance squared
}
// if (Dist2 > 0.01){
for (d = 0; d < dim; ++d){		// Calculate the force
force[d][n] += 24 * (-2*pow(Dist2,-7) + pow(Dist2,-4)) * dist[d];
} 
if (Dist2 < 0.4)
cout << n << " " << m << " " << Dist2 << " " << force[0][1] << endl;
} } }
}

// ________ Update boundaries ________ \\

// copy the particles in the box all around the middle box

void Update_boundaries(){
int b = 0;
for (int x = -1; x < 2; ++x)
for (int y = -1; y < 2; ++y)
for (int z = -1; z < 2; ++z){
for (int n = 0; n < N; ++n){
posb[0][n+(b * N)] = pos[0][n] + boxlength * x; 
posb[1][n+(b * N)] = pos[1][n] + boxlength * y; 
posb[2][n+(b * N)] = pos[2][n] + boxlength * z; 
} 
b++;
}
}


/*
_________________________________________________________

 ------------ Energy related functions ----------------- 
_________________________________________________________

*/

// ________ Calculate Potential Energy ________ \\

// Note: It is more optimal to calculate the force and potential energy in one function

// Calculate potentian energy potential energy only using particles in the box

void Calculate_Potential_Energy(){
double Dist2;
double dist[dim];
EP = 0; 				// Reset total potential energy parameter
for (int n = 0; n < N; ++n){
for (int m = n+1; m < N; ++m){		// Do not count double the potential energy
Dist2 =0;				// Reset distance squared
for (int d = 0; d < dim; ++d){
dist[d] = pos[d][m] - pos[d][n];	// Distance between particles in one 
					// direction
Dist2 += pow(dist[d],2); 		// Distance between particles
}
EP += 4 * (pow(Dist2,-6)-pow(Dist2,-3)); // Total potential energy
} }
}

// Calculate potentian energy potential energy, also using particles in the other boxes

void Calculate_Potential_Energy_with_boundaries(){
double Dist2;
double dist[dim];
EP = 0; 				// Reset total potential energy parameter
for (int n = 0; n < N; ++n){
for (int m = 0; m < 27 * N; ++m){	// Compare with the array containin all parti.
if (13 * N + n != m){			// Only calculate between different 
					// particles (mind that we are 
					// comparing with the boundary matrix 
					// (middle box needed)
Dist2 =0;				// Reset distance squared
for (int d = 0; d < dim; ++d){

dist[d] = posb[d][m] - pos[d][n];	// Distance between particles in one 
					// direction
Dist2 += pow(dist[d],2); 		// Distance between particles
}
EP += 4 * (pow(Dist2,-6)-pow(Dist2,-3))/2; // Total potential energy 
// Note: now I devide it by 2 because the energy will be counted double because it is harder to keep track whether we are dealing with the same particle
} } }
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

 ------------------ Gas quantities --------------------- 
_________________________________________________________

*/

// ________ Calculate the temperature ________ \\

void Calculate_Temperature(){
Temp = (2 * EK) / (kb * N * 3);
}

// ________ Change the temperature ________ \\

// Adjust the velocities of the particles such that the system gains the desired temperature.
// Note: The units of the system are computations units and not natural units

void Adjust_Temperature(){
Calculate_Temperature();	// Use the current temperature
double Temp_Goal = 1;
double avg_vel;
avg_vel = pow((2 * EK) / N,0.5);
for (int n = 0; n < N; ++n){
for (int i = 0; i < dim; ++i){		
vel[i][n] = (Temp_Goal * vel[i][n]) / avg_vel;
} }
}

/*
_________________________________________________________

 ----------------- Extra functions --------------------- 
_________________________________________________________

*/

// ________ Cleaning memory ________ \\

void Clean_up(){
for (int d=0; d<3; ++d){
delete [] pos[d];
delete [] vel[d];
delete [] force[d];
delete [] posb[d];
}
}

// ________ Display values on screen ________ \\

void Display(){
int Display_round = 1;			// Display after "x"  rounds
if (runtemp < 0){
runtemp = Display_round;		// Allways include the first simulation
}
runtemp += 1;
if (runtemp >= Display_round){		// only display every "x" runs
runtemp += -Display_round;  		// reset counter
for (int n = 0; n < N; ++n){
// cout << n << "." << run << " F: " << force[0][n] << " V: " << vel[0][n] << " X: " << pos[0][n] << endl;
// cout << n << "." << run << " X: " << pos[0][n] << " " << "Y: " << pos[1][n] << " " << "Z: " << pos[2][n] << " " << endl;
// cout << n << "." << run << " Fx: " << force[0][n] << " " << "Fy: " << force[1][n] << " " << "Fz: " << force[2][n] << " " << endl;
} 
for (int n = 0; n < 27 * N; ++n){
// cout << n << "." << run << " X: " << posb[0][n] << " " << "Y: " << posb[1][n] << " " << "Z: " << posb[2][n] << " " << endl;
}
cout << "EK: " << EK << " EP: " << EP << " EP + EK: " << EP+EK << endl;
// cout << Temp << endl;
} 
}

// ________ Check for matrix for rows with the same values ________ \\


void Matrix_row_checker(){
for (int n = 0; n < 27 * N; ++n)
for (int m = n + 1; m < 27 * N; ++m)
if (posb[0][n] == posb[0][m])
if (posb[1][n] == posb[1][m])
if (posb[2][n] == posb[2][m])
cout << n << " " << m << endl;
}

/*
_________________________________________________________

 -------------------- Main file ------------------------ 
_________________________________________________________

*/

int main(int argc, char*argv[]){
//--- Read and make variables ---
Parameters(argc, argv);			// Insert parameters
Make_array(argc, argv);			// preallocate arrays
//--- Initiate ---
Initiate_Position_FCC();		// Initiate positions
Initiate_Velocity(); 			// Initiate velocity
//--- Calculate values with initial conditions ---
Update_boundaries();			// surround the box with replica boxes
Calculate_Kinetic_Energy();		// Calculate the kinetic energy
Calculate_Potential_Energy_with_boundaries(); // Calculate the potential energy
		// Alternative: Calculate_Potential_Energy();
for (run = 0; run<Run; ++run){		
//--- Display values ---
Display();				// Display parameters
//--- Update system ---
Update_boundaries();			// surround the box with replica boxes
Update_force_with_boundaries_brute_force(); // Update the force between particles
		// Alternative: Update_force_brute_force();
Update_velocity();			// Update the velocity
Update_position_with_boundaries();	// Update the position
		// Alternative:	Update_position();
//--- Calculate values ---
Calculate_Kinetic_Energy();		// Calculate the kinetic energy
			
//Calculate_Potential_Energy_with_boundaries();
	// Alternative: Calculate_Potential_Energy();
Calculate_Potential_Energy();
Calculate_Temperature();		// Calculate the temperature

// Adjust_Temperature();		// Change the temperature of the system
}


Clean_up();
}
