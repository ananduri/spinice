#include "spinice.h"

void evolve(bool* state, double* intmat, double temp, CRandomMersenne& RanGen_mersenne)
{
	double frac, energy;

	int tau = 5;

	for(int i=0;i<tau;i++)	
	{
		frac = mcstep(RanGen_mersenne, intmat, temp, state, &energy);
		printf("%.3f flipped\n",frac);
	}

	return;
}
