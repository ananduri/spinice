#include "spinice.h"

void evolve(bool* state, int N, double* intmat, double temp, CRandomMersenne& RanGen_mersenne)
{
	double frac;

	int tau;
	if(temp < 1.5)
	{
		tau = 10000;
	}
	else if((temp >= 1.5) && (temp < 2.0))
	{
		tau = 1000;
	}
	else
	{
		tau = 100;
	}

	for(int i=0;i<tau;i++)	
	{
		frac = mcstep(N, RanGen_mersenne, intmat, temp, state);
		printf("%.3f flipped\n",frac);
	}

	return;
}
