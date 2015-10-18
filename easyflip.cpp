#include "spinice.h"

bool easyflip(int N, CRandomMersenne& RanGen_mersenne, bool* spinstate, double* intmat, double T)
{
	extern inline double approxexp(double);

	double r,prob;
	bool flipped = false;	

	int ind = RanGen_mersenne.IRandom(0,N-1);

	double Ediff = 0;
	for(int i=0; i<N; i++)
	{
		Ediff += (spinstate[i]) ? intmat[ind*N + i] : -intmat[ind*N + i];
		//Ediff += (spinstate[i]*2.0 - 1) * intmat[ind*N + i];
	}
	//and what would be faster here? arithmetic, or a ternary statement? i feel like the ternary should be faster,
	//since otherwise it's double arithmetic

	Ediff += (spinstate[ind]) ? -intmat[ind*N + ind] : intmat[ind*N + ind];
	Ediff *= -4;
	
	if(Ediff < 0)
	{
		flipped = true;
		spinstate[ind] ^= true; //flip spin
	}
	else
	{
		r = RanGen_mersenne.Random();
		
		prob = approxexp(-Ediff/T);
		
		if(r < prob)
		{
			flipped = true;
			spinstate[ind] ^= true; //flip spin
		}
	}

	return flipped;
}
