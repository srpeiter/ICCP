#include "Rendering.h"
#include <fstream>
#include <iostream>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

int N;



int main(void)
{
	srand (time(NULL));
	N=40000;
	std::cout << N << ";sjfjsfjs;ljfsl ";
	std::vector<prop> pos(N);



	std::ifstream input("test.txt");


	int n=0;
	while(n < N)
	{

	    input >> pos[n].mass >> pos[n].x >>  pos[n].y >> pos[n].z;

	 n++;
	 
	}
	
	pos[0].xcol = (float)rand()/RAND_MAX;
	std:: cout << rand() << " ";

				    pos[0].ycol = (float)rand()/RAND_MAX;
				    pos[0].zcol = (float)rand()/RAND_MAX;
	for (int n=1 ; n < 10000; n++)
		{pos[n].xcol = pos[n-1].xcol;
			pos[n].ycol =  pos[n-1].ycol;
			    pos[n].zcol =  pos[n-1].zcol;;
		}


	pos[10000].xcol = (float)rand()/RAND_MAX;

				    pos[10000].ycol = (float)rand()/RAND_MAX;
				    pos[10000].zcol= (float) rand()/RAND_MAX;
	for (int n=10001 ; n < 20000; n++)
		{pos[n].xcol = pos[n-1].xcol;
			pos[n].ycol =  pos[n-1].ycol;
			    pos[n].zcol =  pos[n-1].zcol;;
		}

	pos[20000].xcol = (float)rand()/RAND_MAX;

					    pos[20000].ycol = (float)rand()/RAND_MAX;
					    pos[20000].zcol = (float)rand()/RAND_MAX;
		for (int n=20001 ; n < 30000; n++)
			{pos[n].xcol = pos[n-1].xcol;
				pos[n].ycol =  pos[n-1].ycol;
				    pos[n].zcol =  pos[n-1].zcol;;
			}

		pos[30000].xcol = (float)rand()/RAND_MAX;

						    pos[30000].ycol = (float) rand()/RAND_MAX;
						    pos[30000].zcol = (float) rand()/RAND_MAX;
			for (int n=30001 ; n < 40000; n++)
				{pos[n].xcol = pos[n-1].xcol;
					pos[n].ycol =  pos[n-1].ycol;
					    pos[n].zcol =  pos[n-1].zcol;;
				}


	Pre_Render();
	
	int j = 0;
	while(j < 10000)
		Render(pos, N);
	//SaveToTGA();}

	Post_Render();
	
	
	
}
