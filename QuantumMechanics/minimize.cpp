#include"allheaders.h"


// To integrate this code in the program, we should replace the function "dummy" by the energy calculating function


double* minimize(observable& obj,double x0, double x1, double stop_prec)
{
static double output[2];
double eps = x1/100.0;
double x_old = x0;
double f_old = obj.comp_intergral(x0);
double x_new = x1;
double f_new = obj.comp_intergral(x1);
double dx = x_new - x_old;
double df = f_new -f_old;
//To understand the algorithm, let's see what happens in the first iteration
	while(abs(dx) > stop_prec){
	x_old = x_new; // x1 = x_new
    	x_new = x_new - eps * df/dx; // x2 is computed
	dx = x_new - x_old; // dx = x2 -x1 
	f_old = f_new; // f1 = f_new
	f_new = obj.comp_intergral(x_new); // f2 is computed
	df = f_new-f_old; // df = f2-f1
	}

output[0] = x_new;
output[1] = f_new;

return output;
}
