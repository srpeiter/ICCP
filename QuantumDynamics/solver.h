#ifndef SOLVER_H
#define SOLVER_H

#include "math_module.h"
#include<sys/stat.h>

#define pi 3.1415926
#define GNUPLOT "gnuplot -persist "


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

 void setup_grid();

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



#endif
