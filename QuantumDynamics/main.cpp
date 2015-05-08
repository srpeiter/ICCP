#include "solver.h"


#include <unistd.h>
#include <fstream>

#include <signal.h>
#include<sys/stat.h>


int main(void){


double wavevec_x,wavevec_y;
int dimension;
double timestep;
double dist_step ;

std::cout << "Give Dimension of problem" << std::endl;
std:: cout << " Dim=  " ;
std::cin >> dimension;






pid_t gnupid = -3;

//setting gnu plot options

gnupid = fork();

if(gnupid == 0) //child
{

	if (dimension == 1)
	{

std::ofstream cmdfile("gnusettings.txt");

cmdfile   << "set term x11" << "\n"
		  << "set title \" time dependent Schrodinger wavefunction\" " << "\n"
		  << "set time" << "\n"
		  << "set xtics" <<"\n"
		  << "set ytics" << "\n"
		  << "set xrange [0:60]" << "\n"
		  << "set yrange [0:1]" << "\n"
		  //<< "set autoscale y" << "\n"
		  <<"set grid" << "\n"
	      << "plot \'< cat " << "plot.dat" << "\'" << "using 1:2 "<< "with lines"<< "\n"
		  << "pause 0.01 " << "\n"
		  << "reread;" ;
    cmdfile.flush();
	cmdfile.close();
	execlp("gnuplot","gnuplot","gnusettings.txt", (char*) NULL);

	}

	if (dimension == 2)
	{
		std::ofstream cmdfile("gnusettings2d.txt");

		cmdfile   << "set term x11" << "\n"
				  << "set title \" time dependent Schrodinger wavefunction\" " << "\n"
				  << "set time" << "\n"
				  << "set xtics" <<"\n"
				  << "set ytics" << "\n"
				  << "set xrange [0:60]" << "\n"
				  << "set yrange [0:60]" << "\n"
				  << "plot \"plot.dat\" using 1:2:3 with image" << "\n"
				  << "pause 0.01 " << "\n"
				  << "reread;" ;

		    cmdfile.flush();
			cmdfile.close();

			execlp("gnuplot","gnuplot","gnusettings2d.txt", (char*) NULL);
	}




//fork gnuplot process such that it can run independently


}

else if ( gnupid == -1)
		fprintf(stdout," Error forking process");

else{ //parent

mkfifo("plot.dat", S_IWUSR | S_IRUSR); // making a FIFO file to write data for gnuplot

if (dimension == 1)

{

 timestep =0.025;
 dist_step = 0.06;

std::cout << " Give value of wavevector k" << std::endl;
std::cout << " k =  " ;
std:: cin >> wavevec_x;

std::cout << "SIMULATION STARTED" << std::endl;

solver1D Solver(timestep, dist_step, wavevec_x);



Solver.allocate_mem();

Solver.setup_potential();

Solver.init_cond();

Solver.setup_LHS();

Solver.setup_RHS();



Solver.superlu_solve_routine();

Solver.free_mem();

}


if (dimension == 2)
{
	 timestep =0.1;
	 dist_step = 0.3;


	std::cout << " Give values of wavevector k" << std::endl;
	std::cout << " kx =  " ;
	std:: cin >> wavevec_x;

	std::cout << " ky =  " ;
	std:: cin >> wavevec_y;

	std::cout << "SIMULATION STARTED" << std::endl;

	solver2D Solver(timestep, dist_step, wavevec_x ,wavevec_y);


		Solver.allocate_mem();

		Solver.setup_potential();

		Solver.init_cond();



		Solver.setup_lapack_vector_LHS();

		Solver.setup_RHS();


		Solver.lapack_solve_routine();

		Solver.free_mem();


}

kill(gnupid,SIGTERM);

std::cout << "SIMULATION FINISHED" << std::endl;

}

return 0;

}




