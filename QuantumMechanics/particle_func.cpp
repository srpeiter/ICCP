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

// initialize first elements
//

//R1 = {-s/2.0, 0.0, 0.0};
//R2 = { s/2.0, 0.0, 0.0};

R1[0]= -s/2;
R2[0]= s/2;


for (int j=0; j < 3; j++)
{
r1R1_vec[j]= r1[j] - R1[j];
r1R2_vec[j]= r1[j] - R2[j];
r2R1_vec[j]= r2[j] - R1[j];
r2R2_vec[j]= r2[j] - R2[j];
R1R2_vec[j]= R1[j] - R2[j];
r_12vec[j] = r1[j] - r2[j] ;
}

r1R1=0 ; r1R2=0; r2R1=0; r2R2=0; R1R2=0; r_12=0;
//calculating distances;
for (int j=0 ; j < 3 ; j++)
{
r1R1 += r1R1_vec[j]*r1R1_vec[j];
r1R2 += r1R2_vec[j]*r1R2_vec[j];
r2R1 += r2R1_vec[j]*r2R1_vec[j];
r2R2 += r2R2_vec[j]*r2R2_vec[j];
R1R2 += R1R2_vec[j]*R1R2_vec[j];
r_12 += r_12vec[j]*r_12vec[j];

}

}

 
void particle::phi1()
{
phi_1= exp(-r1R1) + exp(-r1R2);
}


void particle::phi2()
{
phi_2= exp(-r2R1) + exp(-r2R2);

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


