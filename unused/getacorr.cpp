#include "spinice.h"

void getacorr(double* acorr, int N, bool* spinstate, double* intmat, double T, int S, int k)
{
	double frac;
	
	double energy;

	srand(time(NULL) + (double)k);
	int seed = rand();
	CRandomMersenne RanGen_mersenne(seed);

	bool* initstate = (bool*)malloc(N*sizeof(bool));

	//equilibrate first
	for(int i=0;i<N;i++)
	{
		initstate[i] = spinstate[i];
	}
	while(calcacorr(N, initstate, spinstate) > 0.01) 
	{
		printf("acorr: %.2f\t",calcacorr(N,initstate,spinstate));
		
		frac=mcstep(N, RanGen_mersenne, intmat, T, spinstate, &energy, k);
		
		printf("equilibrating, %.5f frac flipped\n",frac);
	}


	printf("\n");
	//now get autocorr
	for(int i=0;i<N;i++)
	{
		initstate[i] = spinstate[i];
	}

	acorr[0] = calcacorr(N, initstate, spinstate);
	for(int i=1;i<S;i++)
	{
		printf("acorr: %.3f\tenergy: %.3f\t", acorr[i-1],energy);		

		frac=mcstep(N, RanGen_mersenne, intmat, T, spinstate, &energy, k);

		acorr[i] = calcacorr(N, initstate, spinstate);

		printf("running, %.5f frac flipped\n",frac);
		
	}
	printf("acorr: %.3f\tenergy: %.3f\n",acorr[S-1],energy);

	//repeat, see if acorr is truly the same	
	/*for(int i=0;i<N;i++)
	{
		initstate[i] = spinstate[i];
	}

	acorr[S+0] = calcacorr(N, initstate, spinstate);
	for(int i=1;i<S;i++)
	{
		printf("\tacorr: %.3f\tenergy: %.3f\t", acorr[S+i-1],energy);		

		frac=mcstep(N, RanGen_mersenne, intmat, T, spinstate, &energy, k);

		acorr[S+i] = calcacorr(N, initstate, spinstate);

		printf("\trunning, %.5f frac flipped\n",frac);
		
	}
	printf("\tacorr: %.3f\tenergy: %.3f\n",acorr[S+S-1],energy);

	//repeat again
	for(int i=0;i<N;i++)
	{
		initstate[i] = spinstate[i];
	}

	acorr[S+S+0] = calcacorr(N, initstate, spinstate);
	for(int i=1;i<S;i++)
	{
		printf("\t\tacorr: %.3f\tenergy: %.3f\t", acorr[S+S+i-1],energy);		

		frac=mcstep(N, RanGen_mersenne, intmat, T, spinstate, &energy, k);

		acorr[S+S+i] = calcacorr(N, initstate, spinstate);

		printf("\t\trunning, %.5f frac flipped\n",frac);
		
	}
	printf("\t\tacorr: %.3f\tenergy: %.3f\n",acorr[S+S+S-1],energy);*/


	free(initstate);
	return;
}	
