#include"allheaders.h"
#include<fstream>
#include<iostream>
// this whole code is based on fixing s and varying beta to find minimum.


double observable::energy()
{// here we calculate the enery of the particle
const int dim= 3;
double first_term, second_term, third_term, fourth_term, fifth_term, sixth_term;
double aux1=0, aux2=0, aux3=0, aux4=0, aux5v=0, aux6v=0, aux7v=0, aux8v=0;
double gamma;

double energy;

temp.wavefunction();

gamma = 1 + temp.beta*temp.r_12;



aux1 = (1.0/temp.a) * (exp(-temp.r1R1/temp.a)/temp.phi_1) - 1;
aux2 = (1.0/temp.a) * (exp(-temp.r1R2/temp.a)/temp.phi_1) - 1;
aux3 = (1.0/temp.a) * (exp(-temp.r2R1/temp.a)/temp.phi_2) - 1;
aux4 = (1.0/temp.a) * (exp(-temp.r2R2/temp.a)/temp.phi_2) - 1;


for (int j=0; j < dim; j++)
{
aux5v += (temp.r1R1_vec[j] * temp.r_12vec[j] ) / (temp.r_12*temp.alpha*gamma*gamma);
aux6v += (temp.r1R2_vec[j] * temp.r_12vec[j] ) / (temp.r_12*temp.alpha*gamma*gamma);
aux7v += (temp.r2R1_vec[j] * temp.r_12vec[j] ) / (temp.r_12*temp.alpha*gamma*gamma);
aux8v += (temp.r2R2_vec[j] * temp.r_12vec[j] ) / (temp.r_12*temp.alpha*gamma*gamma);
}



first_term = 1.0/(temp.a*temp.a);

second_term = (1/temp.r_12) * (1 - (2*temp.alpha*gamma + temp.r_12)/(temp.alpha*temp.alpha*pow(gamma,4)));
	
	
third_term = (1/temp.r1R1)* (1 + aux5v)*aux1;
	
	
fourth_term = (1/temp.r1R2)* (1 + aux6v)*aux2;

fifth_term = (1/temp.r2R1) * (1 - aux7v)*aux3;


sixth_term = (1/temp.r2R2) * (1 - aux8v)*aux4;



energy= -first_term + second_term + third_term + fourth_term + fifth_term +sixth_term;

return energy;

}





double observable::comp_integral(double inp_beta, double inp_s)
{
double wave_sq, epsil, fin_energy=0, av_energy;
temp.s=inp_s;
temp.beta=inp_beta;
double stepsize=0.2;
int relax=20000;

int i=0;
temp.generate_metropolis(stepsize);


std:: ifstream infile("metropol.dat");

while (infile >> temp.r1[0] >> temp.r1[1] >> temp.r1[2] >> temp.r2[0] >> temp.r2[1] >> temp.r2[2])
{


if (i > relax)
{
epsil= energy();
fin_energy += epsil;
}

i++;
}

return av_energy = (1.0/(temp.N -relax))*fin_energy;
}


