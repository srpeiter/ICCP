#include"allheaders.h"



double* minimize(particle& obj,double beta0, double beta1, double stop_prec,double s_inp)
{
static double output[2];
double eps = beta1/1000.0;
double beta_old = beta0;

double E_old = obj.comp_integral(beta0,s_inp);

double beta_new = beta1;
double E_new = obj.comp_integral(beta1,s_inp);
double dbeta = beta_new - beta_old;
double dE = E_new -E_old;
// //To understand the algorithm, let's see what happens in the first iteration
// 	while(std::abs(dE) > stop_prec){
// 	beta_old = beta_new; // beta1 = beta_new
//     	beta_new = beta_new - eps * dE/dbeta; // beta2 is computed
// 	dbeta = beta_new - beta_old; // dbeta = beta2 -beta1 
// 	E_old = E_new; // E1 = E_new

// 	E_new = obj.comp_integral(beta_new,s_inp); // E2 is computed
// 	dE = E_new-E_old; // dE = E2-E1
// 	fprintf(stdout,"Beta is: %f,\t Energy is: %f\n", beta_new, E_new);
// 	}

// output[0] = beta_new;
// output[1] = E_new;

double s=1.5;

	for(double beta=0;beta<1;beta=beta+0.03){

		//fprintf(out,"%f\t%f\n",beta,obj.comp_integral(beta,s));
		fprintf(stdout,"%f\t%f\n",beta,obj.comp_integral(beta,s));

	}
return output;
}
