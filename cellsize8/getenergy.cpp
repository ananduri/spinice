#include "spinice.h"

double getenergy(double* intmat, bool* spinstate, double* inter, double* spinstated)
{
	//everytime call this guy, have to create+allocate spinstatei and copy over all the memory
	for(int i=0;i<N;i++)
	{
		spinstated[i] = (spinstate[i]) ? 1.0 : -1.0; 
	}

	cblas_dgemv(CblasColMajor, CblasNoTrans, N, N, 1.0, intmat, N, spinstated, 1, 0.0, inter, 1);

	double sum2 = cblas_ddot(N, spinstated, 1, inter, 1);

	//to keep using the bool array instead, but no BLAS
	/*double sum;
	for(int i=0;i<N;i++)
	{
		sum = 0;
		for(int j=0;j<N;j++)
		{
			sum += intmat[i*N + j] * (2.0*spinstate[j] - 1.0);
		}
		inter[i] = sum;
	}

	double sum2 = 0;
	for(int i=0;i<N;i++)
	{
		sum2 +=	inter[i] * (2.0*spinstate[i] - 1.0);
	}*/

	return sum2;
}
