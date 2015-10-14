#include "spinice.h"

double mcstep(CRandomMersenne& RanGen_mersenne, double* intmat, double T, bool* spinstate, double* energy) {
	
	double counter=0;
	int ind;
	int flip;

	//double* spinstated = (double*)malloc(N*sizeof(double));
	double spinstated[N];

	//double* inter = (double*)malloc(N*sizeof(double)); //allocate this guy once, outside this function?
	double inter[N];

	*energy = getenergy(intmat, spinstate, inter, spinstated);

	for(int i=0;i<N;i++) 
	{
		ind = RanGen_mersenne.IRandom(0,N-1); //inclusive
		
		flip = flipsingle(RanGen_mersenne,ind,spinstate,intmat,T,energy, inter, spinstated);
		
		counter += flip;
	}

	//free(inter);
	//free(spinstated);

	return counter/N;
}
