#include "spinice.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

using namespace std;

int main(int argc, char *argv[]) 
{
	double start,end;
	start = omp_get_wtime();

	//int cellsize = atoi(argv[1]);
	int S = atoi(argv[2]); //number of mc steps
	int step = atoi(argv[3]); 
	double T = strtod(argv[4],NULL);	
	int label = atoi(argv[5]);

	//int N = 16*cellsize*cellsize*cellsize;

	double* intmat = (double*)malloc(N*N*sizeof(double));

	//load intmat
	FILE *matstream;
	char matname[50];
	sprintf(matname,"/home/ananduri/kmc/Intmat_a%d_r3_k10.bin",cellsize);

	matstream = fopen(matname, "rb");
	fread(intmat, sizeof(double), N*N, matstream);
	fclose(matstream);

	srand(time(NULL));
	int seed = rand();
	CRandomMersenne RanGen_mersenne(seed);

	bool* timestate = (bool*)calloc(S*N,sizeof(bool));


	if(step == 0) //start run 
	{
		for(int j=0;j<N;j++)
		{
			timestate[0*N + j] = (RanGen_mersenne.Random() > 0.5);
		}
	}
	else//or read in from file
	{
		FILE* istream;
		char iname[100];
		sprintf(iname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step-1);
	
		istream = fopen(iname,"rb");	
		fseek(istream, -N*sizeof(bool), SEEK_END); 

		fread(timestate, sizeof(bool), N, istream);
		fclose(istream);
	}

	FILE *tstream;
	char tname[100];
	
	//save each one in here
	if(step == 0)  //T and S irrelevant for step = 0
	{
		double Tlim = 2.0;
		double dT = 0.5;
		for (double temp=10;temp>=Tlim;temp-=dT)
		{
			evolve(timestate,intmat,temp,RanGen_mersenne,label);
		}	

		//save
		sprintf(tname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,Tlim,cellsize,label,step); 
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}
	else
	{
		evolvesave(timestate,intmat,T,RanGen_mersenne,S,label);

		//save
		sprintf(tname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step);	
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}


	free(timestate);
	free(intmat);

	printf("cellsize: %d\n\n",cellsize);

	end = omp_get_wtime();
	printf("time: %5.3f\n",end-start);
}
