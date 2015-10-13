#include "spinice.h"

bool flipsingle(CRandomMersenne& RanGen_mersenne, int spinind, bool* spinstate, double* intmat, double T, double* energy, double* inter, double* spinstated)
{
	double r,prob;
	bool flipped = false;	

	double E0 = *energy; //getenergy(N,intmat,spinstate);

	spinstate[spinind] ^= true; //flip spin

	double E1 = getenergy(intmat,spinstate,inter,spinstated);
	
	if( E1 < E0)
	{
		flipped = true;
	}
	else
	{
		r = RanGen_mersenne.Random();
		
		prob = exp(-(E1-E0)/T);

		if(r < prob)
		{
			flipped = true;
		}
		else
		{
			spinstate[spinind] ^= true; //flip back
		}

	}
	
	*energy = (flipped) ? E1 : E0;
	//printf("E1-E0: %.3f\n",E1-E0);
	return flipped;
}
