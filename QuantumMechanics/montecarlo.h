#ifndef MONTECARLO_H
#define MONTECARLO_H

class particle
{

friend double* minimize(particle& obj,double s,int N_beta);

private:		// all variable needed to do computations

double r[6];
double R[6];
double r_12;
double phi_1, phi_2, xii , wave_func;
double alpha, a , beta , s;
int N;

protected:
double r1L_vec[3], r1R_vec[3], r2L_vec[3], r2R_vec[3], r12_vec[3];
double r1L, r1R, r2L, r2R, r12;



public:
particle( double init_pos[6], double beta, double s, int N) :  beta(beta), s(s), N(N)
{ 
  for (int j=0 ; j < 6 ; j++){
				r[j]=init_pos[j]; 
			}
  alpha=2;
}

~particle () {};

void get_a(double s, double criteria);
void initialize();	// initialize everything
void phi1();		
void  phi2();
void xi();
double wavefunction();	// total wavefunction
double energy();	
double comp_integral(double inp_beta, double inp_s);
};

#endif
