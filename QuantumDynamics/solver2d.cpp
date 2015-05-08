#include "solver.h"

void solver2D::writetofile(FILE* pf)
{

	pf=fopen("plot.dat","w");
for (int i = 0; i < N; i++)
	{for(int j = 0 ; j < N ; j++)
	{fprintf(pf, "%f % f  %f  \n",i*dist_step, j*dist_step,  temp_sol[i][j].r*temp_sol[i][j].r + temp_sol[i][j].i*temp_sol[i][j].i);}
	fprintf(pf, "\n");}




	fclose(pf);
}

void solver2D::setup_potential()
{

	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			potential[i][j] = 0 ;
}

void solver2D::init_cond()
{

		double norm=0;
		double rx = 0 , ry = 0 ;

		for(int i = 0 ; i < N; i++)
			for(int j = 0; j < N; j++)
		{
			rx = std::exp(-0.5 * (dist_step*i - (L/2))*(dist_step*i - (L/2)));
			ry = std::exp(-0.5 * (dist_step*j - (L/2))*(dist_step*j - (L/2)));

			wavefunction[i][j].r = rx * ry * std::cos((kx * i * dist_step) + (ky * j * dist_step));	// dont count the boundaries, because they are set to zero
			wavefunction[i][j].i = rx * ry * std::sin((kx * i * dist_step) + (ky * j * dist_step));
		    norm += (( wavefunction[i][j].r*wavefunction[i][j].r + wavefunction[i][j].i*wavefunction[i][j].i)* dist_step*dist_step);
		}

		for(int i= 0; i < N ; i++)
			for(int j = 0 ; j < N; j++)
		{
			wavefunction[i][j].r *= 1/sqrt(norm);
			wavefunction[i][j].i *= 1/sqrt(norm);

		}



}

void solver2D::setup_lapack_vector_LHS()
{
	double g = h_bar*timestep / (4*mass * dist_step*dist_step);
	for(int j= 0 ; j < (N-1) ; j++)
    {
		super_diagonal_LHS[j].i = -g;
        sub_diagonal_LHS[j].i = -g ;
	}

}


// RHS is complex conjugate of LHS
void solver2D::setup_RHS()
{
	double g =  timestep*h_bar/(4*mass*dist_step*dist_step);	//temporary/auxillary variables
											//temporary/auxillary variables

	//making complex matrix

	doublecomplex* row = new doublecomplex [N] ();

	doublecomplex* col= new doublecomplex [N] ();

	// RHS = (1 - i*H)


	row[0].r=1.0 ; row[0].i= -2*g;
	row[1].r=0.0 ; row[1].i= g;

	col[0].r=1.0 ; col[0].i=-2*g;
	col[1].r=0.0 ; col[1].i= g;


math_module::maketoeplitz<doublecomplex>(RHS, row , col ,N);

delete [] row;
delete [] col;
}




void::solver2D::update_RHS(int k, char sweep)
{
	double f = 0;
	if (sweep ==  HOR)
	{


	for(int j = 0 ; j < N ; j++)
	{	f=(timestep*potential[k][j])/(2*pi);
		RHS[j][j].i -= f ;
	}

	}

	else if(sweep == VER)
	{

			for(int j = 0 ; j < N ; j++)
			{	f=(timestep*potential[j][k])/(2*pi);
				RHS[j][j].i -= f ;
			}


	}



}



void solver2D::lapack_update_LHSdiag_RHS(int k , char sweep)
{
	double g = h_bar*timestep / (4*mass * dist_step*dist_step);

	// updating sub and super diagonal of LHS matrix
	for(int j= 0 ; j < (N-1) ; j++)
	    {
			super_diagonal_LHS[j].i = -g;
	        sub_diagonal_LHS[j].i = -g ;
		}


	if ( sweep == HOR)
	{
		update_RHS(k ,sweep);

	for (int j= 0 ; j < N ; j++)
		{
		   diagonal_LHS[j].r = 1.0;
		   diagonal_LHS[j].i = 2*g + (timestep*potential[k][j]/(2*pi));
		}
	}

	else if (sweep == VER)

		update_RHS(k, sweep);

	{
		for (int j= 0 ; j < N ; j++)
				{
				   diagonal_LHS[j].r = 1.0;
				   diagonal_LHS[j].i = 2*g + (timestep*potential[j][k]/(2*pi));
				}



	}


}

void solver2D::lapack_solve_routine()
{
	//fifo file opener
	FILE* pl;

	int info;
	int nrhs = 1;



	for (int k =0 ; k < N ; k++)
	math_module::matvec_multiply<doublecomplex>(out[k], RHS, wavefunction[k], N);

	doublecomplex** temp_sol2 = new doublecomplex*[N];
	for(int i =0 ; i < N; i++)
		temp_sol2[i] = new doublecomplex[N];

//----- loop starts

for ( int j = 0; j < 200; j++) // iteration
{
	//solving horizontal sweep; intermediate solution


	std::cout << " j " << j << std::endl;

	for (int m = 0 ; m < N ; m++)
	{

	 lapack_update_LHSdiag_RHS( m , HOR);

	 zgtsv_(&N, &nrhs, sub_diagonal_LHS, diagonal_LHS, super_diagonal_LHS, out[m], &N, &info );


	}



	math_module::trans_pose<doublecomplex>(temp_sol, out, N);




	for (int k =0 ; k < N ; k++)
		math_module::matvec_multiply<doublecomplex>(temp_sol2[k], RHS, temp_sol[k], N);

	//solving vertical sweep; final solution


		for (int m = 0 ; m < N ; m++)
		{

		lapack_update_LHSdiag_RHS( m , VER);

		zgtsv_(&N, &nrhs, sub_diagonal_LHS, diagonal_LHS, super_diagonal_LHS, temp_sol2[m], &N, &info );

		}

		math_module::trans_pose<doublecomplex>(temp_sol, temp_sol2, N);


	writetofile(pl);





		for (int k =0 ; k < N ; k++)
			math_module::matvec_multiply<doublecomplex>(out [k], RHS, temp_sol[k], N);




	for(int k = 0; k < N; k++)
	{memset(temp_sol[k], 0, N*sizeof(temp_sol[0][0]));
	 memset(temp_sol2[k], 0, N* sizeof(temp_sol[0][0]));
	}


}


//----- loop ends

for (int i = 0 ; i < N; i++)
	delete [] temp_sol2[i];

delete [] temp_sol2;

}
