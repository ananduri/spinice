#include "spinice.h"

bool flipsingle(int N, CRandomMersenne& RanGen_mersenne, int spinind, bool* spinstate, double* intmat, double T, double* energy, double* inter, double* spinstated)
{
	extern inline double approxexp(double);

	double r,prob;
	bool flipped = false;	

	double E0 = *energy; 

	spinstate[spinind] ^= true; //flip spin

	double E1 = getenergy(N,intmat,spinstate,inter,spinstated);
	
	if( E1 < E0)
	{
		flipped = true;
	}
	else
	{
		r = RanGen_mersenne.Random();
		
		prob = approxexp(-(E1-E0)/T);

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
	
	return flipped;
}
