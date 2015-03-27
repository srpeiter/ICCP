#include"allheaders.h"
#include<fstream>
#include<iostream>
// this whole code is based on fixing s and varying beta to find minimum.


double observable::energy()
{// here we calculate the enery of the particle
const int dim= 3;
double energy;

temp.wavefunction();
energy= 0.5*(temp.a*temp.a + temp.x*temp.x*(1- pow(temp.a,4))) ;
	
return energy;

}


double observable::comp_integral(double inp_a){

double epsil, fin_energy=0, av_energy;
temp.a=inp_a;
double stepsize=2;
int relax=100;
int i=0;


srand(time(NULL));

FILE *pf;
double wave_now, wave_next, x_int, condition, rand_num,range;
double x_old;
double delta_E=0;
double E=0;
range=(0.5*RAND_MAX);

pf = fopen("metropol.dat", "w");

for (int i=0 ; i < temp.N ; i++)
{

wave_now=temp.wavefunction();
condition = 0;
rand_num = 0;

rand_num= (double)rand()/(double)RAND_MAX;

x_old=temp.x;



	
		temp.x= x_old+ 2* (((double)rand()/range)-1); 
		wave_next=temp.wavefunction();
		condition=(wave_next*wave_next)/(wave_now*wave_now) ;

	if ( condition < rand_num){
	temp.x= x_old;}
	else{
	delta_E = energy();}

E=E+delta_E;
fprintf(pf,"%f \n",temp.x);
}

fclose(pf);

return E/temp.N;

//std:: ifstream infile("metropol.dat");



//while (infile >> temp.x){

//	if (i > relax){
//		epsil= energy();
//		fin_energy += epsil;
//	}
//	i++;
//}
//return av_energy = (1.0/(temp.N-relax))*fin_energy;
}


