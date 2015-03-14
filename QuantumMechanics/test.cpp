#include "allheaders.h"

using namespace std;

int main(void)

{
 double pos[2][3]= {{0.1,0.23,0.1},{0.5,0.1,0.1}};

particle electron(pos,2.1, 1.2);

//double wave=electron.wavefunction();
observable energy_part(electron, 10000);
//fprintf(stdout, "the wavefunction is %f \n", wave);
double energy=energy_part.energy();

energy_part.norm();
double en= energy_part.comp_integral();
cout << "average energy is "<< en<<endl;
return 0;

}

