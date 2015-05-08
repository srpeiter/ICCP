#include "solver.h"


#include <unistd.h>
#include <fstream>

#include <signal.h>
#include<sys/stat.h>


int main(void){



double kx ,ky ;
kx= -100;
ky= 0;//100;




double timestep =0.1;
double dist_step = 0.3;

pid_t gnupid = -3;

gnupid = fork();

if(gnupid == 0)
{
std::ofstream cmdfile("gnusettings2d.txt");

cmdfile   << "set term x11" << "\n"
		  << "set title \" time dependent Schrodinger wavefunction\" " << "\n"
		  << "set time" << "\n"
		  << "set xtics" <<"\n"
		  << "set ytics" << "\n"
		  << "set xrange [0:60]" << "\n"
		  << "set yrange [0:60]" << "\n"
		  << "set autoscale y" << "\n"

		  << "plot \"plot.dat\" using 1:2:3 with image" << "\n"
		  << "pause 0.01 " << "\n"
		  << "reread;" ;

cmdfile.flush();
cmdfile.close();

//fork gnuplot process such that it can run independently
execlp("gnuplot","gnuplot","gnusettings2d.txt", (char*) NULL);

}

else if ( gnupid == -1)
		fprintf(stdout," Error forking process");

else{



	solver2D Solver(timestep, dist_step, kx ,ky);


	Solver.allocate_mem();

	Solver.setup_potential();

	Solver.init_cond();



	Solver.setup_lapack_vector_LHS();

	Solver.setup_RHS();


	Solver.lapack_solve_routine();

	Solver.free_mem();



    kill(gnupid,SIGTERM);

std::cout << "SIMULATION FINISHED" << std::endl;

}

return 0;

}





