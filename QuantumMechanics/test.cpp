#include "allheaders.h"

using namespace std;

int main(void)

{
 double initposx= 0.3;
double out;
double a=0.1;
const int N=30;
double dat1[N], dat2[N];
double step =0.2;
int it = 100000;
particle electron(initposx, 0.1, it);
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
//out = minimize(energy_part, 2.3, 2.5,0.1,(+0.2*j));
out = energy_part.comp_integral(a+0.1*j);
cout << "a is "<< a<< "and energy is "<<  out << endl;//cout << "average energy is "<< en<<endl;
dat1[j]=(a+0.1*j);
dat2[j]=out;
}

double *out2;
out2=minimize(energy_part, 1.4, 1.8, 0.01 );
cout << " out1 is "<< out2[0] << " and out 2 is "<< out2[1]<< endl;
mydata testclass(dat1, dat2, N);
testclass.plot2d();
return 0;


}

