#include"allheaders.h"




double particle::wavefunction()	// calculating the wavefunction with phi1, phi2 and xi
{
double wave_func= (sqrt(a)/pow(3.14,-4))*exp(-x*x*a*a/2.0);
return wave_func;
}


