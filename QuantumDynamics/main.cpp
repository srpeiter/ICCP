#include "solver.h"


#include <unistd.h>
#include <fstream>



int main(void){

double timestep =0.025;
double dist_step = 0.01;


int gnupid = -3;

//setting gnu plot options

gnupid = fork();

if(gnupid == 0)
{
std::ofstream cmdfile("gnusettings.txt");

cmdfile   << "set title \" time dependent Schrodinger wavefunction\" " << "\n"
		  << "set time" << "\n"
		  << "set xtics" <<"\n"
		  << "set ytics" << "\n"
		  << "set xrange[0:300]" << "\n"
		  << "set yrange[0:0.02]" << "\n"
	      << "plot \'< cat " << "plot.dat" << "\'" <<  "with lines"<< "\n"
		  << "pause 0.1 " << "\n"
		  << "reread;" ;

cmdfile.flush();
cmdfile.close();

//fork gnuplot process such that it can run independently
execlp("gnuplot","gnuplot","gnusettings.txt", (char*) NULL);

}

else if ( gnupid == -1)
		fprintf(stdout," Error forking process");

else{


solver1D Solver(timestep, dist_step);

//Solver.setup_grid();

Solver.allocate_mem();

Solver.setup_potential();

Solver.init_cond();

Solver.setup_RHS();

Solver.setup_LHS();



Solver.superlu_solve_routine();

Solver.free_mem();

}

return 0;

}

/*
const  int N = 5;
const int nnz = 3*N - 2;
doublecomplex** mag = new doublecomplex*[N];
for(int i=0 ; i< N; i++)
	mag[i] = new doublecomplex [N];

doublecomplex* a= new doublecomplex[N] ();
doublecomplex* b= new doublecomplex[N] ();
doublecomplex* data = new doublecomplex[nnz] ();
int* row_ind = new int[nnz] ();
int* col_ptr  = new int[N+1] ();


a[0].r=2.0; a[0].i=3;
a[1].r=-3.0; a[1].i=4;

b[0].r=2.0; b[0].i=3;
b[1].r=-3.0; b[1].i=4;


math_module::maketoeplitz(mag, a,b,N);

SuperMatrix A;

//math_module::CCS_scheme<double>(data, row_ind, col_ptr, mag, N);

math_module::MakeComplexSupermatrix(A,data,row_ind, col_ptr, mag, N);




//math_module::print<double>(mag,N);


/*
std:: cout << "data" << std::endl;

{
	for (int j=0; j < (nnz); j++)
		std::cout << std::setw(3) << data[j].r << "+ " << data[j].i<<' ' ;
}

std:: cout << "row_ind" << std::endl;

{
	for (int j=0; j < (nnz); j++)
		std::cout << std::setw(3) << row_ind[j] << ' ' ;
}

std:: cout << "col_ptr" << std::endl;

{
	for (int j=0; j < (N); j++)
		std::cout << std::setw(3) << col_ptr[j]<< ' ' ;
}




zPrint_CompCol_Matrix("A", &A) ;
delete [] a;
delete [] b;

delete [] col_ptr;
delete [] data;
delete [] row_ind;
for(int i=0 ; i < N; i++)
	delete [] mag[i];

delete [] mag;


return 0;


}

*/


