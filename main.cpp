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
	int step = atoi(argv[3]); 
	double T = strtod(argv[4],NULL);	
	int label = atoi(argv[5]);

	int N = 4*cellsize*cellsize*cellsize;

	double* intmat = (double*)malloc(N*N*sizeof(double));

	//load intmat
	FILE *matstream;
	char matname[50];
	sprintf(matname,"Intmat_a%d_r3_k10.bin",cellsize);

	matstream = fopen(matname, "rb");
	fread(intmat, sizeof(double), N*N, matstream);
	fclose(matstream);

	srand(time(NULL) + label);
	int seed = rand();
	CRandomMersenne RanGen_mersenne(seed);

	bool* timestate = (bool*)calloc(S*N,sizeof(bool));
	
	//LOAD DATA
	if(step == 0 && T >= 1.0) //changed to 1.0, it is now the reference point
	{
		for(int j=0;j<N;j++)
		{
			timestate[0*N + j] = (RanGen_mersenne.Random() > 0.5);
		}
	}
	else if(step == 0 && T < 1.0) //changed to 1.0
	{
		FILE* istream;
		char iname[100];
		sprintf(iname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T+0.1,cellsize,label,step); //now using 0.1 increment, hardcoded 
	
		istream = fopen(iname,"rb");	
		fseek(istream, -N*sizeof(bool), SEEK_END); 

		fread(timestate, sizeof(bool), N, istream);
		fclose(istream);
	}
	else //if step >= 1
	{
		FILE* istream;
		char iname[100];
		sprintf(iname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step-1);
	
		istream = fopen(iname,"rb");	
		fseek(istream, -N*sizeof(bool), SEEK_END); 

		fread(timestate, sizeof(bool), N, istream);
		fclose(istream);
	}


	//EVOLVE DATA
	FILE *tstream;
	char tname[100];

	if(T > 2) //raise error message properly?
	{
		printf("Make T less than 2.0\n");
		exit(1);
	}
	
	double dT;
	if(step == 0 && T >= 1.0) //changed to 1.0, assuming T !> 2
	{
		dT = 0.25;
		for (double temp=10;temp>=2;temp-=dT)
		{
			evolve(timestate,N,intmat,temp,RanGen_mersenne);
		}	

		dT = 0.1;
		for(double temp=2;temp>=T;temp-=dT)
		{
			evolve(timestate,N,intmat,temp,RanGen_mersenne);
		}
		
		evolvesave(timestate,N,intmat,T,RanGen_mersenne,S);

		//save
		sprintf(tname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step); 
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}
	else if(step == 0 && T < 1.0) 
	{
		evolve(timestate,N,intmat,T,RanGen_mersenne);

		evolvesave(timestate,N,intmat,T,RanGen_mersenne,S);	
	
		//save
		sprintf(tname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step);	
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}
	else //if step >=1
	{
		evolvesave(timestate,N,intmat,T,RanGen_mersenne,S);

		//save
		sprintf(tname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step);	
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}


	free(timestate);
	free(intmat);

	printf("cellsize: %d\n\n",cellsize);
	printf("T: %.2f\n\n",T);
	printf("step: %d\n\n",step);

	end = omp_get_wtime();
	printf("wall time: %5.3f\n",end-start);
	printf("in hours: %5.3f\n",(end-start)/3600);
}
