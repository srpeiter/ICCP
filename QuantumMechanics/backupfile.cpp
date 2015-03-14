#include"plotting.h"
#include<random>
#include<chrono>

using namespace std;
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);
uniform_real_distribution<double> distribution(0.0,1.0);

class particle
{

friend class observable; 
 // defining a friend class so to have acces to class particles

private:		// all variable needed to do computations

double r1[3];
double r2[3];
double phi_1L, phi_1R, phi_2L, phi_2R, r_12;
double phi_1, phi_2, xii , wave_func;
double alpha=2, a=0.7 , beta , s;

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

void initialize();	// initialize everything
void phi1();		

void  phi2();

void xi();

double wavefunction();	// total wavefunction


};

void particle::initialize()	// doing the precomputation for the wavefunction and energy
{
for (int j = 0 ; j < 3 ; j++)
{

r_1Lvec[j]= r1[j] + (s/2.0);
r_2Lvec[j]= r2[j] + (s/2.0);
r_1Rvec[j]= r1[j] - (s/2.0);
r_2Rvec[j]= r2[j] - (s/2.0);

r_12vec[j] = r1[j] - r2[j];
}
  
r_1L= r_1Lvec[0]*r_1Lvec[0] + r_1Lvec[1]*r_1Lvec[1] +  r_1Lvec[2]*r_1Lvec[2];
r_2L= r_2Lvec[0]*r_2Lvec[0] + r_2Lvec[1]*r_2Lvec[1] +  r_2Lvec[2]*r_2Lvec[2];
r_1R= r_1Rvec[0]*r_1Rvec[0] + r_1Rvec[1]*r_1Rvec[1] +  r_1Rvec[2]*r_1Rvec[2];
r_2R= r_2Rvec[0]*r_2Rvec[0] + r_2Rvec[1]*r_2Rvec[1] +  r_2Rvec[2]*r_2Rvec[2];

r_12= r_12vec[0]*r_12vec[0] + r_12vec[1]*r_12vec[1] +  r_12vec[2]*r_12vec[2];

r_12 = sqrt(r_12);

phi_1L=exp(-r_1L/a);
phi_1R=exp(-r_1R/a);

phi_2L=exp(-r_2L/a);
phi_2R=exp(-r_2R/a);



}

 
void particle::phi1()
{
phi_1=phi_1L + phi_1R;
}


void particle::phi2()
{
phi_2=phi_2L + phi_2R;
}

void particle::xi()
{
 xii=exp(r_12/ (alpha*(1+beta*r_12)));
}

double particle::wavefunction()	// calculating the wavefunction with phi1, phi2 and xi
{
initialize();
phi1();
phi2();
xi();
return wave_func= phi_1*phi_2*xii;
}


class observable //this class is used for calculating energy integral with montecarlo integration
// This class is a friend of class particle
{
private:
particle temp;
double norm_wave=0;
int N;

public:
observable(particle& temp, int iterations) : temp(temp), N(iterations) {};

double energy();	

void norm();

double comp_integral();

void metropolis_walker();

};

double observable::energy()
{// here we calculate the enery of the particle
const int dim= 3;
double first_term, second_term, third_term, fourth_term, fifth_term, sixth_term;
double inter_res1, inter_res2;
double energy;

temp.wavefunction();

first_term = 1.0/(temp.a*temp.a);
second_term = ((temp.phi_1L/temp.r_1L) + (temp.phi_1R/temp.r_1R))* (1.0/temp.a*temp.phi_1);
third_term =  ((temp.phi_2L/temp.r_2L) + (temp.phi_2R/temp.r_2R))* (1.0/temp.a*temp.phi_2);
fourth_term = (1.0/temp.r_1L) +  (1.0/temp.r_1R) +  (1.0/temp.r_2L) +  (1.0/temp.r_2R);
fifth_term = 1.0/temp.r_12;

for(int j=0; j < dim; j++)
{
inter_res1=((temp.phi_1L*temp.r_1Lvec[j] + temp.phi_1*temp.r_1Rvec[j])/temp.phi_1) - 
((temp.phi_2L*temp.r_2Lvec[j] + temp.phi_2*temp.r_2Rvec[j])/temp.phi_2) ;

inter_res2= temp.r_12vec[j]/(2*temp.a*(1+temp.beta*temp.r_12)*(1+temp.beta*temp.r_12));

fifth_term += inter_res1*inter_res2;
}

sixth_term= ((4*temp.beta +1)*temp.r_12 +4)/(4.0*pow((1+temp.beta*temp.r_12),4)*temp.r_12);

energy= -first_term + second_term + third_term - fourth_term + fifth_term -sixth_term ; 

return energy;

}

void observable::norm()
{//here we normalize the wavefunction only once per s and beta
double wave_sq=0, temp_wave;

for (int j=0; j < N ; j++){

metropolis_walker();
temp_wave= temp.wavefunction();
wave_sq += temp_wave * temp_wave;
}
norm_wave = (1.0/N)* wave_sq;
fprintf(stdout, "the norm of the wave_func is %f\n",norm_wave);

}

double observable::comp_integral()
{
double omega, wave_sq, epsil, fin_energy, av_energy;

norm();
for (int j=0; j < N ; j++)
{

metropolis_walker();

wave_sq= temp.wavefunction();
omega= wave_sq*wave_sq/(norm_wave);
epsil= energy();

fin_energy += omega*epsil;
}
return av_energy = (1.0/N)*fin_energy;
}


void observable::metropolis_walker() // this is the implementation of the a montecarlo method : the metropolis walker
{
double wave_now, wave_next, condition, rand_num;
wave_now=temp.wavefunction();
condition = 0;

temp.r1[0]=distribution(generator);
temp.r1[1]=distribution(generator);
temp.r1[2]=distribution(generator);
temp.r2[0]=distribution(generator);
temp.r2[1]=distribution(generator);
temp.r2[2]=distribution(generator);

wave_next=temp.wavefunction();

condition= (wave_next*wave_next)/(wave_now*wave_now);

if (condition <= 1)
do {
rand_num= distribution(generator);
temp.r1[0]=distribution(generator);
temp.r1[1]=distribution(generator);
temp.r1[2]=distribution(generator);
temp.r2[0]=distribution(generator);
temp.r2[1]=distribution(generator);
temp.r2[2]=distribution(generator);

wave_next=temp.wavefunction();

condition= (wave_next*wave_next)/(wave_now*wave_now);



} while (condition < rand_num);

}


int main(void)

{
 double pos[2][3]= {{0.1,0.23,0.1},{0.5,0.1,0.1}};

particle electron(pos,0.9, 1.2);

//double wave=electron.wavefunction();
observable energy_part(electron, 10000);
//fprintf(stdout, "the wavefunction is %f \n", wave);
double energy=energy_part.energy();

energy_part.norm();
double en= energy_part.comp_integral();
cout << "average energy is "<< en<<endl;
return 0;

}


//double s=1.2;
//double beta= 0.3;
//double r1[3]= {0.3, 0.2, 0.5};
//double r2[3]= {0.5 ,0.1 , 0.1};


