#include"allheaders.h"
#include<fstream>
#include<iostream>
// this whole code is based on fixing s and varying beta to find minimum.


double observable::energy()
{// here we calculate the enery of the particle
const int dim= 3;
double energy;

temp.wavefunction();



energy= temp.a*temp.a + temp.x*temp.x*(1+ pow(temp.a,4)) ;
	
return energy;

}






double observable::comp_integral(double inp_a)
{
double epsil, fin_energy=0, av_energy;
temp.a=inp_a;

double stepsize=0.1;
int relax=100;
int i=0;

temp.generate_metropolis(stepsize);


std:: ifstream infile("metropol.dat");



while (infile >> temp.x)
{
//std::cout << "z " << temp.r2[2] << std::endl;

if (i > relax)
{
epsil= energy();

fin_energy += epsil;
}
i++;
}
return av_energy = (1.0/(temp.N-relax))*fin_energy;
}


