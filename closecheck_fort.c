/* Reads in the output from OUT3 nbody6.
Computes conserved quantities at beginning
and end of run. 
N_tot = first[1], Time = rec[0], N_pairs = rec[1] */

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define N_obj 1001 /* Inital number of objects */
#define N_out 600 // Number of outputs = TCRIT/(DELTAT*NFIX)
#define Obj 922 // Object to be followed
#define G 1.0 /* Newton's constant for AU,yr/2pi,Msun */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define Tconv 0.000697515
#define Rconv 0.000622286
#define Vconv 0.892148


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
        int i,j,k,min_name;
	int N_esc,N_pairs,N_tot;
	FILE *infile;
	int first[7];
	float rec[20];
	float *m,*pos,*vel;
	int *name;
	vector *r,*v,robj,rbh,vobj,vbh;
	float d,min,dcent,vcent,vdiff;

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

		for (j=0;j<N_tot;j++)
		  {
		    r[j].x = pos[3*j];
		    r[j].y = pos[3*j+1];
		    r[j].z = pos[3*j+2];
		    v[j].x = vel[3*j];
		    v[j].y = vel[3*j+1];
		    v[j].z = vel[3*j+2];

		    if(name[j] == Obj) {
			robj.x = r[j].x;
			robj.y = r[j].y;
			robj.z = r[j].z;
			vobj.x = v[j].x;
			vobj.y = v[j].y;
			vobj.z = v[j].z;
		    }
		    if(name[j] == 1) {
			rbh.x = r[j].x;
			rbh.y = r[j].y;
			rbh.z = r[j].z;
			vbh.x = v[j].x;
			vbh.y = v[j].y;
			vbh.z = v[j].z;
		    }
		  }
		dcent = magnitude(vecsub(rbh,robj));
		vcent = magnitude(vecsub(vbh,vobj));
		min = 1000;
		for (j=0;j<N_tot;j++)
		  {
		    // Only particles, not pairs allowed, also not Obj
		      if(name[j] <= N_obj && name[j] != Obj) {
			d = magnitude(vecsub(r[j],robj));
		      }
		      if(d < min) {
			min = d; 
			min_name = name[j];
			vdiff = magnitude(vecsub(v[j],vobj));
		      }
		  }
		rec[0] = rec[0]/(2.0*PI*Tconv);
		min = min/Rconv;
		dcent = dcent/Rconv;
		vcent = vcent/Vconv;
		vdiff = vdiff/Vconv;
		printf("%f %d %f %f %f %f\n",rec[0],min_name,min,dcent,vcent,vdiff);
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
