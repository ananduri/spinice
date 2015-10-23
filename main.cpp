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
	double Jdis = strtod(argv[5],NULL);
	double Ddis = strtod(argv[6],NULL);
	int label = atoi(argv[7]);

	int real_cut,recip_cut;
	if(cellsize==4)
	{
		real_cut = 10;
		recip_cut = 6;
	}
	else if(cellsize==6)
	{
		real_cut = 12;
		recip_cut = 4;
	}
	else
	{
		real_cut =3;
		recip_cut = 10;
	}

	int N = 4*4*cellsize*cellsize*cellsize;

	double* intmat = (double*)malloc(N*N*sizeof(double));

	//load intmat
	FILE *matstream;
	char matname[50];
	if((Jdis != 0.0) && (Ddis != 0.0))
	{
		sprintf(matname,"IntMatdis_J%.2f_D%.2f_a%d_r%d_k%d.bin",Jdis,Ddis,cellsize,real_cut,recip_cut);
	}
	else
	{
		sprintf(matname,"IntMat_a%d_r%d_k%d.bin",cellsize,real_cut,recip_cut);
	}

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
	else if(step == 0 && T < 1.0) 
	{
		FILE* istream;
		char iname[100];
		//sprintf(iname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T+0.1,cellsize,label,step); //now using 0.1 increment, hardcoded 
		sprintf(iname,"samples/a%d/spinJ%.2fD%.2f_T%.2f_a%d_lab%d_step%d.bin",Jdis,Ddis,cellsize,T+0.1,cellsize,label,step); //now using 0.1 increment, hardcoded 
	
		istream = fopen(iname,"rb");	
		fseek(istream, -N*sizeof(bool), SEEK_END); 

		fread(timestate, sizeof(bool), N, istream);
		fclose(istream);
	}
	else //if step >= 1
	{
		FILE* istream;
		char iname[100];
		//sprintf(iname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step-1);
		sprintf(iname,"samples/a%d/spinJ%.2fD%.2f_T%.2f_a%d_lab%d_step%d.bin",Jdis,Ddis,cellsize,T,cellsize,label,step-1); //now using 0.1 increment, hardcoded 
	
		istream = fopen(iname,"rb");	
		fseek(istream, -N*sizeof(bool), SEEK_END); 

		fread(timestate, sizeof(bool), N, istream);
		fclose(istream);
	}


	//EVOLVE DATA
	FILE *tstream;
	char tname[100];
	double dT;

	if(T > 2) //separate for quick high T runs
	{
		dT = 0.25;
		for(double temp=1;temp>=2;temp-=dT)
		{
			evolve(timestate,N,intmat,temp,RanGen_mersenne);
		}
	
		evolvesave(timestate,N,intmat,T,RanGen_mersenne,S);

		//save
		//sprintf(tname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step); 
		sprintf(tname,"samples/a%d/spinJ%.2fD%.2f_T%.2f_a%d_lab%d_step%d.bin",Jdis,Ddis,cellsize,T,cellsize,label,step); 
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}
	
	if(step == 0 && T >= 1.0)
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
		//sprintf(tname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step); 
		sprintf(tname,"samples/a%d/spinJ%.2fD%.2f_T%.2f_a%d_lab%d_step%d.bin",Jdis,Ddis,cellsize,T,cellsize,label,step); 
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}
	else if(step == 0 && T < 1.0) 
	{
		evolve(timestate,N,intmat,T,RanGen_mersenne);

		evolvesave(timestate,N,intmat,T,RanGen_mersenne,S);	
	
		//save
		//sprintf(tname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step);	
		sprintf(tname,"samples/a%d/spinJ%.2fD%.2f_T%.2f_a%d_lab%d_step%d.bin",Jdis,Ddis,cellsize,T,cellsize,label,step); 
		tstream = fopen(tname,"wb");
		fwrite(timestate,sizeof(bool),S*N,tstream);
		fclose(tstream);
	}
	else //if step >=1
	{
		evolvesave(timestate,N,intmat,T,RanGen_mersenne,S);

		//save
		//sprintf(tname,"samples/a%d/spin_T%.2f_a%d_lab%d_step%d.bin",cellsize,T,cellsize,label,step);	
		sprintf(tname,"samples/a%d/spinJ%.2fD%.2f_T%.2f_a%d_lab%d_step%d.bin",Jdis,Ddis,cellsize,T,cellsize,label,step); 
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
