/* Reads in the output from OUT3 nbody6 
Time = rec[7], N pairs = rec[8] */

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define N_obj 101 // Number of objects
#define N_out 500 // Number of outputs = TCRIT/(DELTAT*NFIX)

int main()
{
        int counter,i,N_tot;
	FILE *infile;
	float rec[27];
	float *m,*pos,*vel;
	int *name;
	float comx,comy,comz,mtot = 0;
	float pairs;

	m = malloc(4);
	pos = malloc(4);
	vel = malloc(4);
	name = malloc(4);

       	infile=fopen("OUT3","rb");
       	if (!infile)
       	{
       		printf("Unable to open file!");
       		return 1;
       	}
	for ( counter=0; counter < N_out; counter++)
	{
	        // Read in particle independent data
		fread(&rec,sizeof(rec),1,infile);
		// Adjust total number of objects by number of pairs
		pairs = (int) (rec[8]);
		N_tot = N_obj + pairs;
		// Adjust size of arrays
		m = (float *) realloc(m,N_tot*sizeof(float));
		name = (int *) realloc(name,N_tot*sizeof(int));
		pos = (float *) realloc(pos,3*N_tot*sizeof(float));
		vel = (float *) realloc(vel,3*N_tot*sizeof(float));
		// Read in particle data
       		fread(m,N_tot*sizeof(float),1,infile);
       		fread(pos,3*N_tot*sizeof(float),1,infile);
       		fread(vel,3*N_tot*sizeof(float),1,infile);
       		fread(name,N_tot*sizeof(int),1,infile);
		fseek(infile,4,SEEK_CUR);
		for (i=0;i<N_tot;i++)
		  {
		    if(name[i] == 1) {
		      printf("%d %lg %lg %lg\n",i,pos[3*i],pos[3*i+1],pos[3*i+2]);
		    }
		    comx += m[i]*pos[3*i];
		    comy += m[i]*pos[3*i+1];
		    comz += m[i]*pos[3*i+2];
		    mtot += m[i];
		  }
		comx = comx/mtot;
		comy = comy/mtot;
		comz = comz/mtot;
		//	printf("%lg %lg %lg\n",comx,comy,comz);
	 }
       	fclose(infile);
       	return 0;
}
