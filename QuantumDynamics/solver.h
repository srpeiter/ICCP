#ifndef SOLVER_H
#define SOLVER_H

#include "math_module.h"
#include <sys/stat.h>

#define pi 3.1415926
#define GNUPLOT "gnuplot -persist "
#define  HOR '\x10'
#define  VER '\x20'

// the 1d problem is solved with SUPERLU library
class solver1D
{
protected:
	int nnz;
	doublecomplex** RHS;
	doublecomplex** LHS;
	doublecomplex* wavefunction;
	doublecomplex* out;
	doublecomplex *data;
	double* grid;
	int* row_ind;
	int* col_ptr ;
	int* perm_c;
	int* perm_r;
	double* potential;
	double** hamiltonian;
	double h_bar=1 , mass=1;
	int N;
	double timestep;
	double dist_step;
	double L = 60 ;//+ 2*dist_step ; //length of 1D grid
	double k; //wavevector

public:
solver1D(double inptimestep, double inpdist_step, double wavevector): timestep(inptimestep) , dist_step(inpdist_step) , k(wavevector)
			{ N = (int) (L/dist_step);// -2 ;
			  nnz = 3*N - 2;}  // also count the boundaries (N+2)


void allocate_mem()
{


data = new doublecomplex[nnz] ();
row_ind = new int[nnz] ();
col_ptr = new int[N+1] ();
perm_c = new int[N] ();
perm_r = new int[N] ();
wavefunction = new doublecomplex[N] ();
potential = new double[N] ();
out = new doublecomplex[N] ();
grid = new double[N] ();


RHS = new doublecomplex*[N] ();
for(int i=0; i < N; i++)
	RHS[i]=new  doublecomplex[N] ();

LHS = new doublecomplex* [N] ();
for(int i=0; i < N; i++)
	LHS[i]=new doublecomplex [N]();

hamiltonian = new double* [N] ();
for(int i=0; i < N; i++)
	hamiltonian[i]=new double [N] ();

}
~ solver1D(){};


 void setup_hamiltonian();

 void setup_potential();

 void setup_LHS();

 void setup_RHS();


 void init_cond();

 void writetofile(FILE *pf);


 void  superlu_solve_routine();




void free_mem()
{
delete [] wavefunction;
delete [] potential;
delete [] out;
delete [] perm_c;
delete [] perm_r;
delete [] data;
delete [] row_ind;
delete [] col_ptr;

for(int j= 0 ; j < N; j++)
{
	delete [] RHS[j];
	delete [] LHS[j];
	delete [] hamiltonian[j];
}
delete [] RHS;
delete [] LHS;
delete [] hamiltonian;
delete [] grid;


}

};



// 2d problem is solved with LAPACK library
class solver2D
{
protected:
	int nnz;
	doublecomplex** RHS;

	doublecomplex** temp_sol;
	doublecomplex** wavefunction;
	doublecomplex** out;
	doublecomplex* super_diagonal_LHS;
	doublecomplex* sub_diagonal_LHS;
	doublecomplex* diagonal_LHS;

	double** potential;

	double h_bar=1 , mass=1;
	int N;
	double timestep;
	double dist_step;
	double L = 60 ;//+ 2*dist_step ; //length of 1D grid
	double kx , ky; //wavevectors

public:
solver2D(double inptimestep, double inpdist_step, double wavevector_x, double wavevector_y): timestep(inptimestep) , dist_step(inpdist_step) , kx(wavevector_x) ,
			ky(wavevector_y)
			{ N = (int) (L/dist_step) ;}  // also count the boundaries (N+2)

void allocate_mem()
{


sub_diagonal_LHS = new doublecomplex[N-1] ();
super_diagonal_LHS = new doublecomplex[N-1] ();
diagonal_LHS = new doublecomplex[N] ();


RHS = new doublecomplex*[N] ();
for(int i=0; i < N; i++)
	RHS[i]=new  doublecomplex[N] ();

out = new doublecomplex* [N] ();
for(int i=0; i < N; i++)
	out[i]=new doublecomplex [N]();

temp_sol = new doublecomplex* [N] ();
for(int i=0; i < N; i++)
	temp_sol[i]=new doublecomplex [N] ();

wavefunction = new doublecomplex* [N] ();
for(int i=0; i < N; i++)
	wavefunction[i]=new doublecomplex [N] ();

potential = new double* [N] ();
for(int i=0; i < N; i++)
	potential[i]=new double [N] ();

}

~ solver2D(){};
 void writetofile( FILE* pf);


 void setup_potential();

 void init_cond();

 void setup_lapack_vector_LHS();

 void setup_RHS();


 void update_RHS(int k, char sweep);

 void lapack_update_LHSdiag_RHS(int k, char sweep);



 void  lapack_solve_routine();

 void free_mem()
 {
 delete [] sub_diagonal_LHS;
 delete [] super_diagonal_LHS;
 delete [] diagonal_LHS;


 for(int j= 0 ; j < N; j++)
 {
 	delete [] RHS[j];
 	delete [] temp_sol[j];
 	delete [] potential[j];
 	delete [] out[j];
 	delete [] wavefunction[j];


 }
 delete    [] RHS;
  	delete [] temp_sol;
  	delete [] potential;
  	delete [] out;
  	delete [] wavefunction;

}

};


#endif
