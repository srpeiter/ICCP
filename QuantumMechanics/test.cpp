#include "allheaders.h"

using namespace std;

int main(void)

{
double initposx= 0.3;
double out;
double a=0.1;
const int N=30; 
double dat1[N], dat2[N];
int it = 10000;
double step = 0.2;
particle electron(initposx, 0.1, it);
double an=electron.wavefunction();
electron.generate_metropolis(step);
observable energy_part(electron);

for (int j=0 ; j < N; j++){
	out = energy_part.comp_integral(a+0.1*j);
	cout << "a is "<< a<< " and energy is "<<  out << endl;//cout << "average energy is "<< en<<endl;
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

