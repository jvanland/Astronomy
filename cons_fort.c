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
vector vecadd(vector a, vector b);   /* Vector addition a+b. */
vector vecsub(vector a, vector b);   /* Vector subtraction a-b. */
vector vecmult(float d, vector a);  /* Multiplication: da. */
vector vecdiv(vector a, float d);  /* Division: a/d. */
vector crossprod(vector a, vector b); /* Cross product, axb */
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
	vector *r,*v;
	vector angmom,mom;
	float L,kinetic,potential,tot;
	float d,dcent;

	m = calloc(N_obj,sizeof(float));
	pos = calloc(N_obj,sizeof(float));
	vel = calloc(N_obj,sizeof(float));
	name = calloc(N_obj,sizeof(int));
	r = calloc(N_obj,sizeof(vector));
	v = calloc(N_obj,sizeof(vector));

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
		r = (vector *) realloc(r,N_tot*sizeof(vector));
		v = (vector *) realloc(v,N_tot*sizeof(vector));
		// Read in particle data
       		fread(m,N_tot*sizeof(float),1,infile);
       		fread(pos,3*N_tot*sizeof(float),1,infile);
       		fread(vel,3*N_tot*sizeof(float),1,infile);
       		fread(name,N_tot*sizeof(int),1,infile);
		fseek(infile,4,SEEK_CUR);
	       	if(i==0 || i==N_out-1) { // First and last entries
		  zerovec(&mom);
		  zerovec(&angmom);
		  kinetic=0.0;
		  potential=0.0;

		  for (j=0;j<N_tot;j++)
		    {
		      r[j].x = pos[3*j];
		      r[j].y = pos[3*j+1];
		      r[j].z = pos[3*j+2];
		      v[j].x = vel[3*j];
		      v[j].y = vel[3*j+1];
		      v[j].z = vel[3*j+2];
		      if(name[j] == 1) { // SMBH displacement
			dcent = magnitude(r[j]);
		      }
		      if(name[j] <= N_obj) { // Only particles, not pairs allowed
			mom=vecadd(mom,vecmult(m[j],v[j]));
			angmom=vecadd(angmom,vecmult(m[j],crossprod(r[j],v[j])));
			kinetic+=0.5*m[j]*(v[j].x*v[j].x+v[j].y*v[j].y+v[j].z*v[j].z);
		      }
		    }

		  for (j=0; j<N_tot-1; j++)
		    {
		      for (k=j+1; k<N_tot; k++)
			{
			  if(name[j] <= N_obj && name[k] <= N_obj) {
			    d=sqrt((r[k].x-r[j].x)*(r[k].x-r[j].x)+(r[k].y-r[j].y)*(r[k].y-r[j].y)+(r[k].z-r[j].z)*(r[k].z-r[j].z));
			    potential-=G*m[j]*m[k]/d;
			  }
			}
		    }
		  L=magnitude(angmom);
		  tot = potential+kinetic;
		  printf("%d %f %f %f %f %f %f %f %f\n",N_esc,dcent,tot,potential,kinetic,L,mom.z,mom.y,mom.z);
		}
	 }
       	fclose(infile);
	free(m);
	free(pos);
	free(vel);
	free(name);
	free(r);
	free(v);
       	return 0;
}

void zerovec(vector *r)
{
   (*r).x=0.0;
   (*r).y=0.0;
   (*r).z=0.0;
}

vector vecadd(vector a, vector b)
{
   vector sum;

   sum.x=a.x+b.x;
   sum.y=a.y+b.y;
   sum.z=a.z+b.z;
   return sum;
}

vector vecsub(vector a, vector b)
{
   vector sub;

   sub.x=a.x-b.x;
   sub.y=a.y-b.y;
   sub.z=a.z-b.z;
   return sub;
}

vector vecmult(float d, vector a)
{
   vector mult;

   mult.x=d*a.x;
   mult.y=d*a.y;
   mult.z=d*a.z;
   return mult;
}

vector vecdiv(vector a, float d)
{
   vector div;

   div.x=a.x/d;
   div.y=a.y/d;
   div.z=a.z/d;
   return div;
}

vector crossprod(vector a, vector b)
{
   vector c;

   c.x=a.y*b.z-a.z*b.y;
   c.y=a.z*b.x-a.x*b.z;
   c.z=a.x*b.y-a.y*b.x;
   return c;
}

float magnitude(vector a)
{
   float mag=0.0;

   mag=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
   return mag;
}
