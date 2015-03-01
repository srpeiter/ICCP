#include<iostream>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<random>
#include<unistd.h>

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

char Save_Dir_Energy[] = "Energy.txt";
char Save_Dir_Quantities[] = "Quantities.txt";
double Temp_Goal;
double cutoff = 3;

double EK;		// Total kinetic energy
double EP;		// Total potential energy
double boxlength; 	// The length of the box
double dt;		// discrete time step
double density;		// The denisty of the box
double a;	   	// Length of unit FCC unit cell
double Temp;		// Temperature
double Pressure;	// Pressure
double Temp_old;
/*
double kb = 1.38066*pow(10,-23);// Boltzmann constant
double sigma = 3.405;	// distance to inter-particle distance of zero in natural units (Angstrom)
double epsilon = 119.8 *kb;	// depth of potential well in natural units "kb*T"(K)
double tau = 1.5 * pow(10,-13); // natural time unit
double mass = epsilon * pow(sigma,-2);
*/
double kb = 1;
double sigma = 1;
double epsilon = 1;
double tau = 1;
double mass = 1;



int N;			// Number of particles
int cub_num;		// Number of unit FCC boxes
int Run; 		// Number of simulations rounds
int run;		// Round number
int runtemp = -1;	// Used for the display function 
int bp;			// Number of boundary particles
const int dim = 3;	// Number of dimension

double *pos[dim];	// Position array X,Y,Z
double *posb[dim];	// Position array X,Y,X with boundaries included
			// The size of this matrix is 27 * N
double *vel[dim];	// Velocity array Vx,Vy,Vz
double *force[dim];	// Force array Fx,Fy,Fz
double *force_prev[dim];// Force of the last round


/*
_________________________________________________________

 -------------- Function initialization ---------------- 
_________________________________________________________

*/

// Initiation Functions
void Make_array();					// In use
void Parameters();					// In use		
void Initiate_Position_FCC();				// In use
void Initiate_Velocity();				// In use
void Initiate_Force();					// In use

// Updating functions
void Update_position_with_boundaries_Plus();		// In use 
void Update_velocity_Plus();				// In use 		 
void Update_force_wb_bf_fast();				// In use 	
void Update_boundaries_fast(); 				// In use	

// Energy related functions
void Calculate_Kinetic_Energy_N();			// In use

// Gas quantities
void Calculate_Temperature();				// In use		
void Adjust_Temperature();				// In use

// Combined functions			
void Calculate_PE_Pressure_fast();			// In use	

// Extra functions
void Clean_up();					// In use
void Display();						// In use		
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
force_prev[i]= new double[N] ();
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
//dt = atof(argv[4]);
dt = 0.001;
Temp_Goal = atof(argv[4]);

// Others and calculated values
boxlength = pow((N/density),0.3333); 	// Lengths of the box
cub_num = ceil(pow((N/4),0.333));	// Number of boxes
a=boxlength/cub_num;			// Length of unity box
}
}

// ________ Initiate particles position ________ \\

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

// ________ Initiate particles "Previous" force ________ \\

void Initiate_Force(){
for (int n=0; n<N; ++n){
force[0][n] = 0;
force[1][n] = 0;
force[2][n] = 0;
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

void Update_position_with_boundaries_Plus(){
for (int n = 0; n < N; ++n){
for (int d = 0; d < dim; ++d){		
pos[d][n] += vel[d][n] *dt + force[d][n] * pow(dt,2) /2;	
pos[d][n] += -boxlength * floor(pos[d][n]/boxlength); // subtracts (-1, 0 or 1) * boxlength
}}
}

// ________ Update particles velocity ________ \\

void Update_velocity_Plus(){
for (int n = 0; n < N; ++n){
for (int d = 0; d < dim; ++d){		
vel[d][n] = vel[d][n] + (force[d][n] + force_prev[d][n])*dt/2 ;
}}
}

// ________ Update particles force ________ \\

void Update_force_wb_bf_fast(){
double Dist2;
int d;
double dist[dim]; 
double temp_force;
for (int n = 0; n < N; ++n){
for (d = 0; d < dim; ++d){
force_prev[d][n] = force[d][n];		// save previous force
force[d][n] = 0; 			// Reset the force
}
for (int m = n+1; m < N + bp; ++m)
if (m != n){	// Including boundaries, 27 x N particles
Dist2 = 0;				// Reset "Distance squared"
for (d = 0; d < dim; ++d){
dist[d] = posb[d][m] - pos[d][n];	// Calculate the distance in x,y or z
Dist2 += pow(dist[d],2);		// Calculate the distance squared
}
if ( pow(Dist2,0.5) < cutoff){
for (d = 0; d < dim; ++d){		// Calculate the force
temp_force = 24 * (-2*pow(Dist2,-7) + pow(Dist2,-4)) * dist[d];
force[d][n] += temp_force;
if (m < N){
force[d][m] += -temp_force;
} 
} } } }
}

// ________ Update boundaries ________ \\

void Update_boundaries_fast(){
// Only add particles if neccesary only particles that are withing "cutoff" range of the box
// copy the particles of the box
for (int n = 0; n < N; ++n){
posb[0][n] = pos[0][n]; 
posb[1][n] = pos[1][n]; 
posb[2][n] = pos[2][n]; 
} 
// Add boundary particles
bp = 0;
for (int n = 0; n < N; ++n){
for (int x = -1; x < 2; ++x)
for (int y = -1; y < 2; ++y)
for (int z = -1; z < 2; ++z)
if (x == 0 && y == 0 && z == 0){}
else {
if (	(abs(pos[0][n] + (-0.5 + x)* boxlength) < (0.5 * boxlength + cutoff))&&\
	(abs(pos[1][n] + (-0.5 + y)* boxlength) < (0.5 * boxlength + cutoff))&&\
	(abs(pos[2][n] + (-0.5 + z)* boxlength) < (0.5 * boxlength + cutoff))\
   ){
posb[0][N+bp] = pos[0][n] + boxlength * x; 
posb[1][N+bp] = pos[1][n] + boxlength * y; 
posb[2][N+bp] = pos[2][n] + boxlength * z; 
bp++;
} } }
cout << bp << " " << boxlength << endl;
}

/*
_________________________________________________________

 ------- Energy and Pressure related functions --------- 
_________________________________________________________

*/

// Calculate potentian energy potential energy, also using particles in the other boxes

// ________ Calculate Kinetic Energy ________ \\

void Calculate_Kinetic_Energy(){
EK = 0;      				// Reset total kinetic energy parameter
for (int n = 0; n < N; ++n)
for (int d = 0; d < dim; ++d)
EK += 0.5 * pow(vel[d][n],2);	 // Calculate kinet energy /particle/dim
}

void Calculate_PE_Pressure_fast(){
// Calculate the pressure and potential energy

double Dist2;
double dist[dim];
Pressure = density * Temp;
EP = 0; 				// Reset total potential energy parameter
for (int n = 0; n < N; ++n){
for (int m = 0; m < N + bp; ++m){	// Compare with the array containin all parti.
if (n != m){				// Only calculate between different 
					// particles (mind that we are 
					// comparing with the boundary matrix 
					// (middle box needed)
Dist2 =0;				// Reset distance squared
for (int d = 0; d < dim; ++d){
dist[d] = posb[d][m] - pos[d][n];	// Distance between particles in one 
					// direction
Dist2 += pow(dist[d]/sigma,2); 	// Distance between particles
}
if ( pow(Dist2,0.5) < cutoff){
EP += epsilon * 4 * (pow(Dist2,-6)-pow(Dist2,-3))/2; // Total potential energy 
// Note: now I devide it by 2 because the energy will be counted double because it is harder to keep track whether we are dealing with the same particle
if (m > n){
for (int d = 0; d < dim; ++d){		// Calculate the force
Pressure += 24 * (-2*pow(Dist2,-6) + pow(Dist2,-3)) ;
} } } } } }
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
// Calculate_Temperature();	// Use the current temperature

Temp_old = Temp;

for (int n = 0; n < N; ++n){
for (int i = 0; i < dim; ++i){		
vel[i][n] = vel[i][n] * pow(Temp_Goal / Temp,0.5);
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
// cout << Temp << " " << Pressure << endl;
} 
}


// ________ Print data ________ \\

class graphdata
{
int particles;
double *xdata=NULL;
double *ydata=NULL;
double *zdata=NULL;
public:
graphdata (double *xdata, int particles) : xdata(xdata), particles(particles) {}
graphdata (double *xdata, double *ydata, int particles): xdata(xdata), ydata(ydata), particles(particles) {}
graphdata (double *xdata, double *ydata, double *zdata, int particles): xdata(xdata), ydata(ydata), zdata(zdata), particles(particles) {}
void printtofile(char filename[]);
};
void graphdata::printtofile(char filename[])
{
FILE *f;
f=fopen(filename, "w");
if (xdata != NULL && ydata == NULL)
{
for (int i =0; i < particles; i++)
fprintf(f, "%f\n",xdata[i]);
}
if (xdata != NULL && ydata != NULL && zdata ==NULL)
{
for (int i =0; i < particles; i++)
fprintf(f, "%f %f\n",xdata[i], ydata[i]);
}
if (xdata != NULL && ydata != NULL && zdata !=NULL)
{
for (int i =0; i < particles; i++)
fprintf(f, "%f %f %f\n",xdata[i], ydata[i], zdata[i]);
}
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
double EP_data[Run];
double EK_data[Run];
double Pressure_data[Run];
double Temp_data[Run];
double Temp_dif[Run];

//--- Initiate ---
Initiate_Position_FCC();		// Initiate positions
Initiate_Velocity(); 			// Initiate velocity
Initiate_Force();			// Initiate (indirectly) the previous  
					// forceses by initializing the force as 0
//--- Calculate values with initial conditions ---
Update_boundaries_fast();		// surround the box with replica boxes
Calculate_Kinetic_Energy();		// Calculate the kinetic energy
Calculate_PE_Pressure_fast();
for (run = 0; run<Run; ++run){		
//--- Display values ---
Display();				// Display parameters
//--- Update system ---
Update_force_wb_bf_fast();
Update_velocity_Plus();			// Update the velocity
Update_position_with_boundaries_Plus();	// Update the position
Update_boundaries_fast();		// surround the box with replica boxes
// cout << "bp: " << bp << endl;
//--- Calculate values ---
Calculate_Kinetic_Energy();		// Calculate the kinetic energy		
Calculate_PE_Pressure_fast();
Calculate_Temperature();		// Calculate the temperature

//--- Save data ---
EK_data[run] = EK;
EP_data[run] = EP;
Pressure_data[run] = Pressure;
Temp_data[run] = Temp;
Temp_dif[run] = Temp - Temp_old;
cout << Temp << endl;
//--- Adjust system --- 
Adjust_Temperature();		// Change the temperature of the system

}
//--- Print data ---
graphdata Energy(EP_data, EK_data,Run);
Energy.printtofile(Save_Dir_Energy);
graphdata Quantities(Temp_data,Pressure_data,Temp_dif,Run);
Quantities.printtofile(Save_Dir_Quantities);
//--- Clean memory ---
Clean_up();
}
