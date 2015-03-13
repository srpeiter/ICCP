#include"plotting.h"


using namespace std;

class particle
{

friend class part_prop;

private:
double r1[3];
double r2[3];
double phi_1L, phi_1R, phi_2L, phi_2R, r_12;
double alpha=2, a=1 , beta , s;

protected:
double r_1L, r_1R, r_2L, r_2R;
double r_12vec[3], r_1Lvec[3], r_1Rvec[3], r_2Lvec[3], r_2Rvec[3];


public:
particle( double position[][3], double beta, double s) :  beta(beta), s(s)
{ r1[0] = position[0][0];
  r1[1] = position[0][1];
  r1[2] = position[0][2]; 
  r2[0] = position[1][0];
  r2[1] = position[1][1];
  r2[2] = position[1][2];
  
  pre_comp();
}

~particle () {};

void pre_comp();
double phi1();

double phi2();

double xi();

double wavefunction();


};

void particle::pre_comp()
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

 
double particle::phi1()
{
return phi_1L + phi_1R;
}


double particle::phi2()
{
return phi_2L + phi_2R;
}

double particle::xi()
{
return exp(r_12/ (alpha*(1+beta*r_12)));
}

double particle::wavefunction()
{
return phi1() * phi2() * xi();
}



int main(void)

{
 double pos[2][3]= {{0.1,0.2,0.3},{0.5,0.3,0.2}};

particle electron(pos,0.2, 1.2);

double wave=electron.wavefunction();

fprintf(stdout, "the wavefunction is %f \n", wave);
cout << "wave is "<< wave <<endl;
return 0;

}


