#include "spinice.h"

double mcstep(int N, CRandomMersenne& RanGen_mersenne, double* intmat, double T, bool* spinstate) {
	
	//int ind;

	//double* spinstated = (double*)malloc(N*sizeof(double));

	//double* inter = (double*)malloc(N*sizeof(double)); //allocate this guy once, outside this function?

	//*energy = getenergy(N, intmat, spinstate, inter, spinstated); //initialize energy; from now on only do simple method to get it
	
	int flip;

	double counter=0;
	for(int i=0;i<N;i++) 
	{
		//ind = RanGen_mersenne.IRandom(0,N-1); //inclusive //combine this into easyflip
		
		//flip = flipsingle(N,RanGen_mersenne,ind,spinstate,intmat,T,energy,inter,spinstated);
		flip = easyflip(N,RanGen_mersenne,spinstate,intmat,T);
		
		counter += flip;
	}

	//free(inter);
	//free(spinstated);

	return counter/N;
}
