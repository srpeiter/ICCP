#include "solver.h"


void solver1D:: setup_grid()
{
//for(int i=0 ; i < N; i++)
	//grid[i]= i *dist_step;

}


void solver1D::writetofile(FILE* pf)
{
	pf =fopen("plot.dat","w");
	for (int j=0 ; j < N; j++)
	fprintf(pf, "%f  %f \n", j*dist_step, out[j].r*out[j].r + out[j].i*out[j].i);
	fclose(pf);
}



void solver1D::setup_hamiltonian()
{
// make tridiagonal matrix
double g =  timestep*h_bar/(4*mass*dist_step*dist_step);	//temporary/auxillary variables
double f = 0;											//temporary/auxillary variables

// toeplitz setup arrays
double row[N]={0};
double col[N]= {0};

row[0]= 2*g;
row[1]= -g;

col[0]= 2*g;
col[1]=-g;

math_module::maketoeplitz<double>(hamiltonian, row , col , N);

//include time-independent potential in hamiltonian
for(int i = 0 ; i < N ; i++)
{	f=(timestep*potential[i])/(2*pi);
	hamiltonian[i][i] += f;
}

}

void solver1D::setup_potential()
{
double omega = 100;
for(int i=0 ; i < N ; i++)
	potential[i]=0.5*(i*dist_step - (L/2)) * (i*dist_step - (L/2));// *omega;

}

void solver1D::init_cond()
{
	FILE* pfile;
	pfile = fopen("test.dat","w");
	double norm=0;
	double k = 1;
	for(int j = 0 ; j < int(N); j++)
	{	wavefunction[j].r = std::exp(-0.5 * (dist_step*j - (L/2))*(dist_step*j - (L/2))) * std::cos(k * j * dist_step);	// dont count the boundaries, because they are set to zero
		wavefunction[j].i = std::exp(-0.5 * (dist_step*j - (L/2))*(dist_step*j - (L/2))) * std::sin(k * j * dist_step);
	norm += (( wavefunction[j].r*wavefunction[j].r + wavefunction[j].i*wavefunction[j].i)* j * dist_step);}

	for(int i= 0; i < N ; i++)
	{
		wavefunction[i].r *= 1/sqrt(norm);
		wavefunction[i].i *= 1/sqrt(norm);}



}


void solver1D::setup_LHS()
{

double g =  timestep*h_bar/(4*mass*dist_step*dist_step);	//temporary/auxillary variables
double f = 0;											//temporary/auxillary variables

//making complex matrix

doublecomplex* row = new doublecomplex [N] ();

doublecomplex* col= new doublecomplex [N] ();

// RHS = (1 - i*H)


row[0].r=1.0 ; row[0].i= 2*g;
row[1].r=0.0 ; row[1].i= -g;

col[0].r=1.0 ; col[0].i=2*g;
col[1].r=0.0 ; col[1].i= -g;



math_module::maketoeplitz<doublecomplex>(RHS, row , col ,N);
for(int i = 0 ; i < N ; i++)
{	f=(timestep*potential[i])/(2*pi);
	RHS[i][i].i += f ;
}

delete [] row;
delete [] col;

}

void solver1D::setup_RHS()
{
// LHS complex conjugate of RHS
for(int i=0 ; i < N ; i++)
	for(int j=0; j < N; j++)
	{	LHS[i][j].r = RHS[j][i].r;
		LHS[i][j].i = -RHS[j][i].i;}
}



void solver1D::superlu_solve_routine()
{


// SUPERLU library options
SuperMatrix A;
nnz = 3*N - 2;
int nrhs = 1;
superlu_options_t options;
SuperLUStat_t stat;
int info;
set_default_options(&options);
options.ColPerm = NATURAL;
StatInit(&stat);
//-------------------------

// plotting setup-------------------------

mkfifo("plot.dat", S_IWUSR | S_IRUSR);
FILE* pl;
//---------------------------------

// Solving --------------------------


math_module::MakeComplexSupermatrix<doublecomplex>(A, data, row_ind, col_ptr, LHS, N );

math_module::matvec_multiply<doublecomplex>(out, RHS, wavefunction , N);

doublecomplex* temp_sol = new doublecomplex[N] ();   // temporary buffer


for(int j=0 ; j < 1000; j++)
{

SuperMatrix  B, L, U;
zCreate_Dense_Matrix(&B, N, nrhs, out, N, SLU_DN, SLU_Z, SLU_GE);

// Solve the linear system

zgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);


writetofile(pl); // put data into fifo


math_module::matvec_multiply<doublecomplex>(temp_sol, RHS, out , N);

memcpy(out, temp_sol, N*sizeof(temp_sol[0]));
memset(temp_sol, 0, sizeof(temp_sol[0])*N);

Destroy_SuperMatrix_Store(&B);
Destroy_SuperNode_Matrix(&L);
Destroy_CompCol_Matrix(&U);







}

StatFree(&stat);
delete [] temp_sol;


// -------------------------------




}



