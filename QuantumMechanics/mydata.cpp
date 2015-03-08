
#include<iostream>
#include<stdio.h>
#include<stdlib.h>

using namespace std;
int N=20;
double xmeas[20];
double ymeas[20];
double zmeas[20];


void maketestdata()
{
	for (int i; i<N; i++)
	{
	xmeas[i]=i;
	ymeas[i]=3*i;
	zmeas[i]=4*i+4;

	}



}



class mydata
{
int particles;
double *xdata=NULL;
double *ydata=NULL;
double *zdata=NULL;
public:
mydata (double *xdata, int particles) : xdata(xdata), particles(particles) {}

mydata (double *xdata, double *ydata, int particles): xdata(xdata), ydata(ydata), particles(particles) {}

mydata (double *xdata, double *ydata, double *zdata, int particles): xdata(xdata), ydata(ydata), zdata(zdata), particles(particles) {}


void printtofile(char filename[]);
void plot2d();
void plot3d();
void histplot();
};


void mydata::printtofile(char filename[])
{
FILE *f;
f=fopen(filename, "w");
if (xdata != NULL && ydata == NULL)
{
for (int i =0; i < particles; i++)
fprintf(f, "%f\n",xdata[i]);
}
if (xdata != NULL && ydata != NULL && zdata ==NULL)
{
for (int i =0; i < particles; i++)
fprintf(f, "%f %f\n",xdata[i], ydata[i]);
}
if (xdata != NULL && ydata != NULL && zdata !=NULL)
{
for (int i =0; i < particles; i++)
fprintf(f, "%f %f %f\n",xdata[i], ydata[i], zdata[i]);
}
}


void mydata::plot2d()
{
printtofile("test.dat");
system ("./myplot test.dat using 1:2");

}


int main()
{
maketestdata();

mydata Testdata(xmeas,ymeas,N);

Testdata.plot2d();

return 0;


}
