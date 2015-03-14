#include"allheaders.h"
#include<chrono>

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::uniform_real_distribution<double> distribution(0.0,1.0);




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



