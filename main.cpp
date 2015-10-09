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

	int cellsize = atoi(argv[1]);
	int S = atoi(argv[2]); //number of mc steps
	int step = atoi(argv[3]); //number of separate runs, useless now
	double T = strtod(argv[4],NULL);	
	int label = atoi(argv[5]);

	int N = 16*cellsize*cellsize*cellsize;

	double* intmat = (double*)malloc(N*N*sizeof(double));

	//load intmat
	FILE *matstream;
	char matname[50];
	sprintf(matname,"Intmat_a%d_r3_k10.bin",cellsize);

	matstream = fopen(matname, "rb");
	fread(intmat, sizeof(double), N*N, matstream);
	fclose(matstream);

	srand(time(NULL));
	int seed = rand();
	CRandomMersenne RanGen_mersenne(seed);

	/*double** acorr = (double**)malloc(M*sizeof(double*));
	for(int i=0;i<M;i++) {
		acorr[i] = (double*)calloc(S,sizeof(double)); //delete the 3* later
	}*/


	/*bool** timestate = (bool**)malloc(S*sizeof(bool*));
	for(int i=0;i<S;i++)
	{
		timestate[i] = (bool*)calloc(N,sizeof(bool));
	}*/

	bool* timestate = (bool*)calloc(S*N,sizeof(bool));



	if(step == 0) //start run 
	{
		//intialize with shitty RNG
		for(int j=0;j<N;j++)
		{
			//timestate[0*N + j] = ((double)rand()/(double)RAND_MAX > 0.5) ? true : false;
			//timestate[0*N + j] = ((double)rand()/(double)RAND_MAX > 0.5) ? 1 : -1;
			timestate[0*N + j] = (RanGen_mersenne.Random() > 0.5);
		}
	}
	else//or read in from file
	{
		FILE* istream;
		char iname[100];
		sprintf(iname,"samples/a%d/spin_T%.2f_S%d_a%d_lab%d_step%d.bin",cellsize,T,S,cellsize,label,step-1);
	
		istream = fopen(iname,"rb");	
		fseek(istream, -N*sizeof(bool), SEEK_END); 

		fread(timestate, sizeof(bool), N, istream);
		fclose(istream);
	}
	

	FILE *tstream;
	char tname[100];
	//put in annealing schedule here
	//save each one, in here
	//also need reading them in
	double dT = 0.5;
	for (double temp=10;temp>=2;temp-=dT)
	{
		evolve(timestate,N,intmat,temp,RanGen_mersenne,label);

		evolvesave(timestate,N,intmat,temp,RanGen_mersenne,S,label);

		//save
		sprintf(tname,"samples/a%d/spin_T%.2f_S%d_a%d_lab%d_step%d.bin",cellsize,temp,S,cellsize,label,step);	
		//sprintf(tname,"samples/a%d/Ac_T%.2f_S%d_a%d_lab%d_step%d.txt",cellsize,temp,S,cellsize,label,step);	
		
		/*tstream = fopen(tname, "w");
		for(int j=0;j<S;j++)
		{
			for(int m=0;m<N;m++)
			{
				fprintf(tstream,  "%d ", timestate[j][m]);
			}
			fprintf(tstream, "\n");
		}
		fclose(tstream);*/

		//write binary instead, test this
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}

	/*
	//for M trials,
	//#pragma omp parallel num_threads(1) 
	{
	bool* spinstate = (bool*)malloc(N*sizeof(bool));
	//#pragma omp for
	for(int i=0;i<M;i++) 
	{
		// initialize spin array with shitty RNG
		for(int j=0;j<N;j++) {
			spinstate[j] = ((double)rand()/(double)RAND_MAX > 0.5) ? true : false;
		}

		// equilibrate and calculate autocorrelation
		getacorr(acorr[i],N,spinstate,intmat,T,S,(i*10)+label);

		printf("run %d done\n",i);
		printf("\n");
	}
	free(spinstate);
	}*/


	free(intmat);

	//FILE *tstream;
	//char tname[50];
	//sprintf(tname,"tstate_




	/*for(int i=0;i<S;i++)
	{
		free(timestate[i]);
	}*/
	free(timestate);

	

	//save all data, for each trial
	/*FILE *acstream;
	char mcname[50];
	sprintf(mcname,"Ac_T%.2f_S%d_a%d_lab%d.txt",T,S,cellsize,label);	

	acstream = fopen(mcname, "w");
	for(int i=0;i<M;i++) {
		for(int j=0;j<S;j++) { //also a 3* here
			fprintf(acstream, "%.5f ",acorr[i][j]);
		}
		fprintf(acstream,"\n");
	}
	fclose(acstream);

	for(int i=0;i<M;i++) {
		free(acorr[i]);
	}
	free(acorr);*/

	//printf("cellsize: %d\nT: %.2f\n\n",cellsize,T);
	printf("cellsize: %d\n\n",cellsize);

	end = omp_get_wtime();
	//double time = ((double)(end-start))/CLOCKS_PER_SEC;
	printf("time: %5.3f\n",end-start);
}
