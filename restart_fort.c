/* Reads in the output from OUT3 nbody6.
Creates a new start file from an intermediate time point.
N_tot = first[1], Time = rec[0], N_pairs = rec[1] */

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

int main(int argc, char *argv[ ])
{
  if(argc != 5) {
    printf("Usage: %s N_obj N_out N_res Obj\n", argv[0]);
    return 1;
  }
        int i,j;
	int N_tot;
	int N_obj = atoi(argv[1]);
	int N_out = atoi(argv[2]);
	int Restart = atoi(argv[3]);
	int Obj = atoi(argv[4]);
	FILE *infile,*outfile;
	int first[7];
	float rec[20];
	float *m,*pos,*vel;
	int *name;

	m = malloc(4);
	pos = malloc(4);
	vel = malloc(4);
	name = malloc(4);

       	infile=fopen("OUT3","rb");
	outfile=fopen("fort.restart","w");
       	if (!infile)
       	{
       		printf("Unable to open file!");
       		return 1;
       	}
	for (j=0;j<N_out;j++)
	{
	        // Read in particle independent data
		fread(&first,sizeof(first),1,infile);
		fread(&rec,sizeof(rec),1,infile);
		// Determine N_tot
		N_tot = first[1];
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
		if (j == Restart) {
		  for (i=0;i<N_tot;i++) {
		    if (name[i] <= N_obj) {  // Pairs not allowed
		      if(name[i] == Obj) printf("Obj %d at position %d\n",Obj,i+1);
		      fprintf(outfile,"%f %f %f %f %f %f %f\n",m[i],pos[3*i],pos[3*i+1],pos[3*i+2],vel[3*i],vel[3*i+1],vel[3*i+2]);
		    }
		  }
		}
	 }
       	fclose(infile);
	free(m);
	free(pos);
	free(vel);
	free(name);
       	return 0;
}
