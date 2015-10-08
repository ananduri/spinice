#include "spinice.h"

double mcstep(int N, CRandomMersenne& RanGen_mersenne, double* intmat, double T, bool* spinstate, double* energy, int k) {
	
	double counter=0;
	int ind;
	int flip;

	//srand(time(NULL) + (double)k);
	//int seed = rand();
	//CRandomMersenne RanGen_mersenne2(seed);

	*energy = getenergy(N, intmat, spinstate);

	for(int i=0;i<N;i++) 
	{
		ind = RanGen_mersenne.IRandom(0,N-1); //if use this numbers are same everytime
		//printf("flipind: %d\n",ind); 
		
		flip = flipsingle(N,RanGen_mersenne,ind,spinstate,intmat,T,energy);
		//flip = flipsingle(N,RanGen_mersenne,ind,spinstate,intmat,T);
		
		counter += flip;
	}

	return counter/N;
}
