#include "allheaders.h"

using namespace std;

int main(void)

{
 double pos[2][3]= {{0.1,0.23,0.1},{0.5,0.1,0.1}};
double *out;
double s=0.1;
const int N=30;
double dat1[N], dat2[N];
double step =0.2;
int it = 100000;
particle electron(pos,2.1,0.5,it);
double an=electron.wavefunction();
electron.generate_metropolis(step);
observable energy_part(electron);

//double wave=electron.wavefunction();
//fprintf(stdout, "the wavefunction is %f \n", wave);
//double energy=energy_part.comp_integral(0.4,1);
//fprintf(stdout, "the energy is %f \n", energy);

//energy_part.norm();
//double en= energy_part.comp_integral(0.9,0.5);


for (int j=0 ; j < N; j++)
{
out = minimize(energy_part, 8.8, 9,0.1,(s+0.1*j));
cout << "beta is "<< out[0]<< " and energy is "<<out[1]<< endl;
//cout << "average energy is "<< en<<endl;
dat1[j]=(s+0.1*j);
dat2[j]=out[1];
}

mydata testclass(dat1, dat2, N);
testclass.plot2d();
return 0;


}

