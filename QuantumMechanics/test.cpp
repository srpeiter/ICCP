#include "allheaders.h"

using namespace std;

int main(void)

{
double initposx= 0.1;
double out;
double a_start=0.5;
double a_end=1.5;
double a_step=0.01;
const int N=(a_end-a_start)/a_step; 
double dat1[N], dat2[N];
int it = 100000;


particle electron(initposx, 0.1, it);
double an=electron.wavefunction();
electron.generate_metropolis(2);
observable energy_part(electron);

for (int j=0 ; j < N; j++){
	out = energy_part.comp_integral(a_start+a_step*j);
	cout << "a is "<< a_start+a_step*j<< " and energy is "<<  out << endl;//cout << "average energy is "<< en<<endl;
	dat1[j]=(a_start+a_step*j);
	dat2[j]=out;
}

double *out2;
out2=minimize(energy_part, 1.4, 1.8, 0.01 );
cout << " out1 is "<< out2[0] << " and out 2 is "<< out2[1]<< endl;
mydata testclass(dat1, dat2, N);
testclass.plot2d();
return 0;


}

