
class mydata
{
int particles;
double *xdata=NULL;
double *ydata=NULL;
double *zdata=NULL;
public:
graphdata (double *xdata, int particles) : xdata(xdata), particles(particles) {}
graphdata (double *xdata, double *ydata, int particles): xdata(xdata), ydata(ydata), particles(particles) {}
graphdata (double *xdata, double *ydata, double *zdata, int particles): xdata(xdata), ydata(ydata), zdata(zdata), particles(particles) {}
void printtofile(char filename[]);
};
void graphdata::printtofile(char filename[])
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


