/* Reads in the output from OUT3 nbody6.
Computes conserved quantities at beginning
and end of run. 
N_tot = first[1], Time = rec[0], N_pairs = rec[1] */

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define N_obj 101 /* Inital number of objects */
#define N_out 5000 // Number of outputs = TCRIT/(DELTAT*NFIX)
#define G 1.0 /* Newton's constant for AU,yr/2pi,Msun */
#define M 1.0e2    /* Mass of SMBH in solar masses */
#define PI 3.141592653589793  /* pi to unnecessary precision */


typedef struct
{
   float x;     /* x component of vector. */
   float y;     /* y component of vector. */
   float z;     /* z component of vector. */
} vector;

void zerovec(vector *r);    /* Set components of vector to zero */
vector vecsub(vector a, vector b);   /* Vector subtraction a-b. */
float magnitude(vector a);  /* Magnitude of vector a */

int main()
{
        int i,j,k;
	int N_esc,N_pairs,N_tot;
	FILE *infile;
	int first[7];
	float rec[20];
	float *m,*pos,*vel;
	int *name;
	vector r1,r2;
	float d;

	m = calloc(N_obj,sizeof(float));
	pos = calloc(N_obj,sizeof(float));
	vel = calloc(N_obj,sizeof(float));
	name = calloc(N_obj,sizeof(int));

       	infile=fopen("OUT3","rb");
       	if (!infile)
       	{
       		printf("Unable to open file!");
       		return 1;
       	}
	for (i=0;i<N_out;i++)
	{
	        // Read in particle independent data
		fread(&first,sizeof(first),1,infile);
		fread(&rec,sizeof(rec),1,infile);
		// Determine N_esc,N_pairs,N_tot
		N_tot = first[1];
		N_pairs = (int) rec[1];
		N_esc = N_obj - (N_tot - N_pairs);
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
	       	if(i==N_out-1) { // Last entry
		  for (j=0;j<N_tot;j++)
		    {
		      if(name[j] == 1) {
			r1.x = pos[3*j];
			r1.y = pos[3*j+1];
			r1.z = pos[3*j+2];
		      }
		    }
		  for (j=0;j<N_tot;j++) {
		    r2.x = pos[3*j];
		    r2.y = pos[3*j+1];
		    r2.z = pos[3*j+2];
		    d = magnitude(vecsub(r1,r2));
		      printf("%d %f\n",name[j],d);
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

void zerovec(vector *r)
{
   (*r).x=0.0;
   (*r).y=0.0;
   (*r).z=0.0;
}

vector vecsub(vector a, vector b)
{
   vector sub;

   sub.x=a.x-b.x;
   sub.y=a.y-b.y;
   sub.z=a.z-b.z;
   return sub;
}

float magnitude(vector a)
{
   float mag=0.0;

   mag=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
   return mag;
}
