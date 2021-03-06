#include"plotting.h"

// function of the plotting library 
// plotting done in gnuplot
// see plot2d, plot3d and hist for small scripts written in gnuplot
void mydata::printtofile(const char filename[])
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
system ("./plot2d");

}

void mydata::histplot()
{
printtofile("test.dat");
system("./hist");
}

void mydata::plot3d()
{
printtofile("test.dat");
system("./plot3d");

}

