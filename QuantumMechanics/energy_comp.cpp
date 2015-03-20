#include"allheaders.h"
#include<fstream>
// this whole code is based on fixing s and varying beta to find minimum.


double observable::energy()
{// here we calculate the enery of the particle
const int dim= 3;
double first_term, second_term, third_term, fourth_term, fifth_term, sixth_term, seventh_term;
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
inter_res1=((temp.phi_1L*temp.r_1Lvec[j] + temp.phi_1R*temp.r_1Rvec[j])/temp.phi_1) - 
((temp.phi_2L*temp.r_2Lvec[j] + temp.phi_2R*temp.r_2Rvec[j])/temp.phi_2) ;

inter_res2= temp.r_12vec[j]/(2*temp.a*(1+temp.beta*temp.r_12)*(1+temp.beta*temp.r_12));

sixth_term += inter_res1*inter_res2;
}

seventh_term= ((4*temp.beta +1)*temp.r_12 +4)/(4.0*pow((1+temp.beta*temp.r_12),4)*temp.r_12);

energy= -first_term + second_term + third_term - fourth_term + fifth_term +sixth_term - seventh_term; 

return energy;

}

void observable::norm()
{//here we normalize the wavefunction only once per s and beta
double wave_sq=0, temp_wave;

std:: ifstream infile("metropol.dat");

while (infile >> temp.r1[0] >> temp.r1[1] >> temp.r1[2] >> temp.r2[0] >> temp.r2[1] >> temp.r2[2])
{

temp_wave= temp.wavefunction();
wave_sq += temp_wave * temp_wave;
}
norm_wave = (1.0/N)* wave_sq;
fprintf(stdout, "the norm of the wave_func is %f\n",norm_wave);

}





double observable::comp_integral(double inp_beta, double inp_s)
{
double omega, wave_sq, epsil, fin_energy, av_energy;
temp.s=inp_s;
temp.beta=inp_beta;
double stepsize=0.1;
int relax=10000;
temp.generate_metropolis(stepsize);
norm();

std:: ifstream infile("metropol.dat");

while (infile >> temp.r1[0] >> temp.r1[1] >> temp.r1[2] >> temp.r2[0] >> temp.r2[1] >> temp.r2[2])
{

wave_sq= temp.wavefunction();
omega= wave_sq*wave_sq/(norm_wave);
epsil= energy();


fin_energy += omega*epsil;
}
return av_energy = (1.0/(N))*fin_energy;
}


