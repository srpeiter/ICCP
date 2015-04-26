#include "solver.h"


#include <unistd.h>
#include <fstream>

#include <signal.h>
#include<sys/stat.h>


int main(void){

double wavevector ;

std::cout << "Give a value for the wavevector k" << std::endl;
std:: cout << "k= " << std::endl;

std::cin >> wavevector;

std::cout << "STARTING SIMULATION ....." << std::endl;

double timestep =0.025;
double dist_step = 0.06;


pid_t gnupid = -3;

//setting gnu plot options

gnupid = fork();

if(gnupid == 0)
{
std::ofstream cmdfile("gnusettings.txt");

cmdfile   << "set term x11" << "\n"
		  << "set title \" time dependent Schrodinger wavefunction\" " << "\n"
		  << "set time" << "\n"
		  << "set xtics" <<"\n"
		  << "set ytics" << "\n"
		  << "set xrange [0:60]" << "\n"
		  << "set yrange [0:0.002]" << "\n"
		  //<< "set autoscale y" << "\n"
		  <<"set grid" << "\n"
	      << "plot \'< cat " << "plot.dat" << "\'" << "using 1:2 "<< "with lines"<< "\n"
		  << "pause 0.01 " << "\n"
		  << "reread;" ;

cmdfile.flush();
cmdfile.close();

//fork gnuplot process such that it can run independently
execlp("gnuplot","gnuplot","gnusettings.txt", (char*) NULL);

}

else if ( gnupid == -1)
		fprintf(stdout," Error forking process");

else{


solver1D Solver(timestep, dist_step, wavevector);

Solver.setup_grid();

Solver.allocate_mem();

Solver.setup_potential();

Solver.init_cond();

Solver.setup_LHS();

Solver.setup_RHS();



Solver.superlu_solve_routine();

Solver.free_mem();

kill(gnupid,SIGTERM);

std::cout << "SIMULATION FINISHED" << std::endl;

}

return 0;

}




