#include "spinice.h"

void evolvesave(bool* timestate, int N, double* intmat, double T, CRandomMersenne& RanGen_mersenne, int S, int label)
{
	double frac, energy;

	bool* tempstate = (bool*)malloc(N*sizeof(bool));
	for(int i=0;i<N;i++)
	{
		tempstate[i] = timestate[i];
	}
	
	for(int i=0;i<S-1;i++)
	{
		frac = mcstep(N, RanGen_mersenne, intmat, T, tempstate, &energy, label);
	
		for(int j=0;j<N;j++)
		{
			timestate[(i+1)*S + j] = tempstate[j];
		}
	}

	free(tempstate);

	return;
}
