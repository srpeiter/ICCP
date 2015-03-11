
#include"plotting.h"

using namespace std;
const int N=1000;
double xmeas[N];
double ymeas[N];
double zmeas[N];

default_random_engine generator;
normal_distribution<double> norm_dist(0.0,1.0) ;


void maketestdata()
{
	for (int i; i<N; i++)
	{
	xmeas[i]=norm_dist(generator);
	ymeas[i]=norm_dist(generator);
	zmeas[i]=norm_dist(generator);

	}



}

int main()
{
maketestdata();

mydata Testdata(xmeas,N);
Testdata.printtofile("test.dat");

//Testdata.plot2d();
Testdata.histplot();
//Testdata.plot3d();
return 0;


}
