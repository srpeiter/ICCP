#ifndef MONTECARLO_H
#define MONTECARLO_H

class particle
{

friend class observable;
 // defining a friend class so to have acces to class particles


private:		// all variable needed to do computations

double r1[3],R1[3] ={0};
double r2[3],R2[3] = {0};
double r_12;
double phi_1, phi_2, xii , wave_func;
double alpha, a , beta , s;
int N;

protected:
double r1R1_vec[3], r1R2_vec[3], r2R1_vec[3], r2R2_vec[3], r_12vec[3], R1R2_vec[3];
double r1R1, r1R2, r2R1, r2R2, R1R2;


public:
particle( double position[][3], double beta, double s, int N) :  beta(beta), s(s), N(N)
{ r1[0] = position[0][0];		//initializing positions
  r1[1] = position[0][1];
  r1[2] = position[0][2]; 
  r2[0] = position[1][0];
  r2[1] = position[1][1];
  r2[2] = position[1][2];
  alpha=2;
}

~particle () {};

void get_a(double& s, double criteria);
void initialize();	// initialize everything
void phi1();		

void  phi2();

void xi();

double wavefunction();	// total wavefunction

void generate_metropolis(double step);



};

class observable //this class is used for calculating energy integral with montecarlo integration
// This class is a friend of class particle
{
// minimize is a friend function of class observable
friend double* minimize(observable& obj, double x0, double x1, double stop_prec, double inp_s );

private:
particle temp;

public:
observable(particle& temp) : temp(temp) {}


double energy();	


void metropolis_walker();

double comp_integral(double inp_beta, double inp_s);
};



#endif
