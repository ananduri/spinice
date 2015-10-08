#include "spinice.h"

void evolve(bool* state, int N, double* intmat, double temp, CRandomMersenne& RanGen_mersenne,int label)
{
	double frac, energy;
	//int tau = schedule[(int)temp*2];
	//int tau = schedule[0]; 
	int tau = 20;

	for(int i=0;i<tau;i++)	
	{
		frac = mcstep(N, RanGen_mersenne, intmat, temp, state, &energy, label);
	}

	return;
}
