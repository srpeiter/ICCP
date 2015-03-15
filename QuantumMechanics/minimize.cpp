#include"allheaders.h"


// To integrate this code in the program, we should replace the function "dummy" by the energy calculating function


double* minimize(observable& obj,double x0, double x1, double stop_prec,double s_inp)
{
static double output[2];
double eps = x1/1000.0;
double x_old = x0;

double f_old = obj.comp_integral(x0,s_inp);
fprintf(stdout,"the old energy is %f\n", f_old);

double x_new = x1;
double f_new = obj.comp_integral(x1,s_inp);
double dx = x_new - x_old;
fprintf(stdout,"delta beta is  %f\n", dx);
double df = f_new -f_old;
//To understand the algorithm, let's see what happens in the first iteration
	while(std::abs(dx) > stop_prec){
	x_old = x_new; // x1 = x_new
    	x_new = x_new - eps * df/dx; // x2 is computed
	dx = x_new - x_old; // dx = x2 -x1 
	f_old = f_new; // f1 = f_new
fprintf(stdout,"x_new is  %f\n", x_new);

	f_new = obj.comp_integral(x_new,s_inp); // f2 is computed
	fprintf(stdout,"the new energy is %f\n", f_new);
	df = f_new-f_old; // df = f2-f1
	}

output[0] = x_new;
output[1] = f_new;

return output;
}
