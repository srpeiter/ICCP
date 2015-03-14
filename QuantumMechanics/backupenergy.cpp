#include<random>
#include<math.h>
#include"plotting.h"

using namespace std;
const int N=1000;
double energy[N];

/* This code computes the Energy 
 * by integrating a 6-D integral.
 *  
 */

void compute_integral(double &beta, double &s, double *r1, double *r2)
{
for (int i=0 ; i < N ; i++)
{
const int dim= 3;
double first_term, second_term, third_term, fourth_term, fifth_term, sixth_term;
double a=1.1;
double phi_1, phi_2, phi_1L, phi_1R, phi_2L, phi_2R;
double r_12, r_1L, r_1R, r_2L, r_2R;
double r_12vec[3], r_1Lvec[3], r_1Rvec[3], r_2Lvec[3], r_2Rvec[3];
double inter_res1=0, inter_res2=0, inter_res3=0;

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

phi_1L=exp(-r_1L/a);
phi_1R=exp(-r_1R/a);

phi_2L=exp(-r_2L/a);
phi_2R=exp(-r_2R/a);

phi_1= phi_1L + phi_1R;
phi_2= phi_2L + phi_2R;

	
first_term = 1.0/(a*a);
second_term = ((phi_1L/r_1L) + (phi_1R/r_1R))* (1.0/a*phi_1);
third_term =  ((phi_2L/r_2L) + (phi_2R/r_2R))* (1.0/a*phi_2);
fourth_term = (1.0/r_1L) +  (1.0/r_1R) +  (1.0/r_2L) +  (1.0/r_2R);
fifth_term = 1.0 / abs(r_12);

for(int j=0; j < dim; j++)
{
inter_res1=((phi_1L*r_1Lvec[j] + phi_1*r_1Rvec[j])/phi_1) - 
((phi_2L*r_2Lvec[j] + phi_2*r_2Rvec[j])/phi_2) ;

inter_res2= r_12vec[j]/(2*a*(1+beta*r_12)*(1+beta*r_12));

fifth_term += inter_res1*inter_res2;
}

sixth_term= ((4*beta +1)*r_12 +4)/(4.0*pow((1+beta*r_12),4)*r_12);

energy[i]= -first_term + second_term + third_term - fourth_term + fifth_term -sixth_term ; 

}

}

int main()
{

double s=1.2;
double beta= 0.3;
double r1[3]= {0.3, 0.2, 0.5};
double r2[3]= {0.5 ,0.1 , 0.1};

compute_integral(beta , s, r1, r2);


mydata Testdata(energy,N);

Testdata.printtofile("test.dat");

//Testdata.plot2d();
//Testdata.histplot();
//Testdata.plot3d();
return 0;


}
