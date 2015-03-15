#include"allheaders.h"


void particle::get_a(double& s, double criteria)
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


void particle::initialize()	// doing the precomputation for the wavefunction and energy
{
get_a(s, 0.1);	//compute a first

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


