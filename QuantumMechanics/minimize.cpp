#include"allheaders.h"



double* minimize(particle& obj,double s,int N_beta){

double beta_min = 0.3;
double beta_max = 1;

int N=N_beta;
static double output[2];
double x[N];
double y[N];

double Sx=0,Sx2=0,Sy=0;
double Sxx=0,Sxy=0,Sxx2=0,Sx2y=0,Sx2x2=0;
double out[4];
double a,b,c,r2;
double Sresid=0,Stotal=0,yresid;

do{
for (int i = 0; i < N; ++i)
{
	x[i] = beta_min + (double)(i+1)*(beta_max-beta_min)/(double)N;
	y[i] = obj.comp_integral(x[i],s);
	fprintf(stdout, "", y[i]);
}

for (int i = 1; i < N-2; ++i)
{
	if((y[i]>y[i+2])&&(y[i]>y[i-1])){
		y[i] = (y[i+2]+y[i-1])/2;
	}
}


for (int i = 0; i < N; ++i)
{
	Sx 		+=x[i];
	Sx2 	+=x[i]*x[i];
	Sy 		+=y[i];

	Sxx 	+=x[i]*x[i];
	Sxy 	+=x[i]*y[i];
	Sxx2 	+=x[i]*x[i]*x[i];
	Sx2y 	+=x[i]*x[i]*y[i];
	Sx2x2 	+=x[i]*x[i]*x[i]*x[i];
}

Sxx 	 	-=Sx*Sx/N;
Sxy 		-=Sx*Sy/N;
Sxx2  		-=Sx*Sx2/N;
Sx2y  		-=Sx2*Sy/N;
Sx2x2 		-=Sx2*Sx2/N;

a = (Sx2y*Sxx-Sxy*Sxx2)/(Sxx*Sx2x2-Sxx2*Sxx2);
b = (Sxy*Sx2x2-Sx2y*Sxx2)/(Sxx*Sx2x2-Sxx2*Sxx2);
c = Sy/N - b*Sx/N - a*Sx2/N;

for (int i = 0; i < N; ++i)
{
	yresid = y[i]- (a*x[i]*x[i] + b*x[i] + c);
	Sresid += yresid*yresid;

	Stotal +=(y[i]-Sy/N)*(y[i]-Sy/N)/N;
}

r2 = 1-Sresid/Stotal;

output[0] = a*b/2/a*b/2/a - b*b/2/a + c;
output[1] = - b/2/a;

fprintf(stdout, "energy: %f\n",output );

}while(a<0);


return output;}
