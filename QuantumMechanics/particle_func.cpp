#include"allheaders.h"

// this function calculates "a" for a given "s". 
void particle::get_a(double s, double criteria)
{
	double l = 0.5;
	double r = 1;
	a = 0.75;
	double a_old = 0;

	while (std::abs(a_old-a) > criteria){
		(a*log(a/(1-a))>s) ? (r=a) : (l=a);	
		a_old = a;	
		a = (r+l)/2;
	}
}

// this function is called every time "s" changes 
// It inializes everything and it precomputes the variable needed to calcuate the wavefucntions
void particle::initialize()	// doing the precomputation for the wavefunction and energy
{
get_a(s, 0.1);	//compute a first

// initialize first elements

R[0]= -s/2.0;
R[1]=0.;
R[2]=0.;
R[3]= s/2.0;
R[4]=0.;
R[5]=0.;


for (int j=0; j < 3; j++)
{
	r1L_vec[j]= r[j] 	- R[j];
	r1R_vec[j]= r[j] 	- R[j+3];
	r2L_vec[j]= r[j+3] 	- R[j];
	r2R_vec[j]= r[j+3] 	- R[j+3];
	r12_vec[j] = r[j] 	- r[j+3] ;
}

r1L=0 ; r1R=0; r2L=0; r2R=0; r12=0;
//calculating distances;
for (int j=0 ; j < 3 ; j++)
{
	r1L += r1L_vec[j]*r1L_vec[j];
	r1R += r1R_vec[j]*r1R_vec[j];
	r2L += r2L_vec[j]*r2L_vec[j];
	r2R += r2R_vec[j]*r2R_vec[j];
	r12 += r12_vec[j]*r12_vec[j];
}

r1L = std::sqrt(r1L);
r1R = std::sqrt(r1R);
r2L = std::sqrt(r2L);
r2R = std::sqrt(r2R);
r12 = std::sqrt(r12);

}

//wavefunction of particle 1
void particle::phi1()
{
	
	phi_1= exp(-r1L/a) + exp(-r1R/a);
}

//wavefunction of particle 2
void particle::phi2()
{
	phi_2= exp(-r2L/a) + exp(-r2R/a);

}
// interacting wavefunction between particle 1 and 2
void particle::xi()
{
	xii=exp(r12/ (alpha*(1+beta*r12)));
}

//total wavefunction
double particle::wavefunction()	// calculating the wavefunction with phi1, phi2 and xi
{
	initialize();
	phi1();
	phi2();
	xi();
	wave_func= phi_1*phi_2*xii;
	return wave_func;
}

// function energy is just a implementation of the local energy, so it is 
// wise to comment the functioning of every variable, just go through the code 
// and try to understand it

double particle::energy()
{// here we calculate the energy of the particle
	double term_1, term_2, term_3, term_4, term_5, term_6, term_7, term_8; 
	double energy=0;
	wavefunction();

	double phi_1L = exp(-r1L/a);
	double phi_1R = exp(-r1R/a);
	double phi_2L = exp(-r2L/a);
	double phi_2R = exp(-r2R/a);

	term_1 = -1/a/a;
	term_2 = 1/a/phi_1*(phi_1L/r1L+phi_1R/r1R);
	term_3 = 1/a/phi_2*(phi_2L/r2L+phi_2R/r2R);
	term_4 = -1/r1L-1/r1R-1/r2L-1/r2R;
	term_5 = 1/r12;

	term_6 = 0;
	for (int j=0 ; j < 3 ; j++){
	term_6 += ((phi_1L*r1L_vec[j]/r1L+phi_1R*r1R_vec[j]/r1R)/phi_1-(phi_2L*r2L_vec[j]/r2L+phi_2R*r2R_vec[j]/r2R)/phi_2)*r12_vec[j]/(r12*2*a*(1+beta*r12)*(1+beta*r12));
	}

	term_7 = -((4*beta+1)*r12+4)/(4*(1+beta*r12)*(1+beta*r12)*(1+beta*r12)*(1+beta*r12)*r12);
	term_8 = 1/s;

	energy = term_8 + term_7 + term_6 + term_5 + term_4 + term_3 + term_2 + term_1;

	return energy;
}

// this is part where the metropolis walker start making its steps
// go sequentially through the code

double particle::comp_integral(double inp_beta, double inp_s)
{

	double wave_now, wave_next, rand_num,range,delta_E,E; // temporary buffers 
	double condition = 0;
	double r_old[6];
	s=inp_s;
	beta=inp_beta;		// parameters we want to minimize for every "s"
	int relax=(int)2*N/10;  // thermalization steps: bring the system to "equillibrium

	// random number generator settings
	srand(time(NULL)); //seed
	range=(0.5*RAND_MAX);

	for (int j=0 ; j < 6 ; j++){
			r[j]=1; 
		}

	for (int i=0 ; i < N ; i++)
	{

		wave_now=wavefunction();
		rand_num= (double)rand()/(double)RAND_MAX;

		for (int j=0 ; j < 6 ; j++){
			r_old[j]=r[j]; 
		}

		for (int j=0 ; j < 6 ; j++){
			r[j]=r_old[j]+ 2*  (((double)rand()/range)-1); 
		}

		wave_next=wavefunction();
		condition=(wave_next*wave_next)/(wave_now*wave_now) ;

		if( condition < rand_num){

			for (int j=0 ; j < 6 ; j++){
				r[j]=r_old[j]; 
			}

		}
		else{
			if (i>relax)  // thermalization
			{
				delta_E = energy();
			}
		}
		E = E + delta_E;
	}

	return E/(N-relax); // local energy
}
