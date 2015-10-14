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

	int S = atoi(argv[1]); //number of mc steps
	int step = atoi(argv[2]); 
	double T = strtod(argv[3],NULL);	
	int label = atoi(argv[4]);

	double* intmat = (double*)malloc(N*N*sizeof(double));

	//load intmat
	FILE *matstream;
	char matname[50];
	sprintf(matname,"/home/ananduri/kmc/Intmat_a%d_r3_k10.bin",cellsize);

	matstream = fopen(matname, "rb");
	fread(intmat, sizeof(double), N*N, matstream);
	fclose(matstream);

	srand(time(NULL) + label);
	int seed = rand();
	CRandomMersenne RanGen_mersenne(seed);

	bool* timestate = (bool*)calloc(S*N,sizeof(bool)); //declare on stack?
	

	if(step == 0 && T >= 2.0) //start run 
	{
		for(int j=0;j<N;j++)
		{
			timestate[0*N + j] = (RanGen_mersenne.Random() > 0.5);
		}
	}
	else if(step == 0 && T < 2.0)
	{
		FILE* istream;
		char iname[100];
		sprintf(iname,"/home/ananduri/kmc/samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T+0.1,cellsize,label,step); //0.1 increment hardcoded in 
	
		istream = fopen(iname,"rb");	
		fseek(istream, -N*sizeof(bool), SEEK_END); 

		fread(timestate, sizeof(bool), N, istream);
		fclose(istream);
	}
	else//or read in from file
	{
		FILE* istream;
		char iname[100];
		sprintf(iname,"/home/ananduri/kmc/samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step-1);
	
		istream = fopen(iname,"rb");	
		fseek(istream, -N*sizeof(bool), SEEK_END); 

		fread(timestate, sizeof(bool), N, istream);
		fclose(istream);
	}

	FILE *tstream;
	char tname[100];
	
	//save each one in here
	if(step == 0 && T >= 2.0)  
	{
		double dT = 0.5;
		for (double temp=10;temp>=T;temp-=dT)
		{
			evolve(timestate,intmat,temp,RanGen_mersenne);
		}	
		
		evolvesave(timestate,intmat,T,RanGen_mersenne,S);

		//save
		sprintf(tname,"/home/ananduri/kmc/samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step); 
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}
	else if(step == 0 && T < 2.0)
	{
		evolvesave(timestate,intmat,T,RanGen_mersenne,S);	
	
		sprintf(tname,"/home/ananduri/kmc/samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step);	
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}
	else
	{
		evolvesave(timestate,intmat,T,RanGen_mersenne,S);

		//save
		sprintf(tname,"/home/ananduri/kmc/samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step);	
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
