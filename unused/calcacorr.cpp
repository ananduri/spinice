#include "spinice.h"

double calcacorr(int N, bool* initstate, bool* spinstate) {

	double dip1[3];	
	double dip2[3];
	
	double ac = 0;
	
	int i;

	for(i=0;i<N;i++)
	{
		dipassign(dip1, i % 4, initstate[i]);
		dipassign(dip2, i % 4, spinstate[i]);
	
		ac += dip1[0]*dip2[0] + dip1[1]*dip2[1] + dip1[2]*dip2[2];
	}

	return ac/N;
}
