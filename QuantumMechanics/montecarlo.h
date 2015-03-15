#ifndef MONTECARLO_H
#define MONTECARLO_H

class particle
{

friend class observable; 
 // defining a friend class so to have acces to class particles


private:		// all variable needed to do computations

double r1[3];
double r2[3];
double phi_1L, phi_1R, phi_2L, phi_2R, r_12;
double phi_1, phi_2, xii , wave_func;
double alpha=2, a , beta , s;

protected:
double r_1L, r_1R, r_2L, r_2R;
double r_12vec[3], r_1Lvec[3], r_1Rvec[3], r_2Lvec[3], r_2Rvec[3];


public:
particle( double position[][3], double beta, double s) :  beta(beta), s(s)
{ r1[0] = position[0][0];		//initializing positions
  r1[1] = position[0][1];
  r1[2] = position[0][2]; 
  r2[0] = position[1][0];
  r2[1] = position[1][1];
  r2[2] = position[1][2];

}

~particle () {};

void get_a(double& s, double criteria);
void initialize();	// initialize everything
void phi1();		

void  phi2();

void xi();

double wavefunction();	// total wavefunction


};

class observable //this class is used for calculating energy integral with montecarlo integration
// This class is a friend of class particle
{
// minimize is a friend function of class observable
friend double* minimize(observable& obj, double x0, double x1, double stop_prec, double inp_s );

private:
particle temp;
double norm_wave=0;
int N;

public:
observable(particle& temp, int iterations) : temp(temp), N(iterations) {};

double energy();	

void norm();

void metropolis_walker();

double comp_integral(double inp_beta, double inp_s);
};



#endif
