/* Converts Cartesian coordinates in bt 
format to Orbital elements. Follows a 
single object. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 39.47841760435743 /* Newton's constant for AU,yr,Msun */
#define M 1.0e5    /* Mass of SMBH in solar masses */
#define N_obj 10000     /* Number of stars */
#define N_steps 18000000    /* Number of pkd steps */
#define Out_int 3000    /* Interval of pkd output in steps */
#define PI 3.141592653589793  /* pi to unnecessary precision */


typedef struct
{
   double x;     /* x component of vector. */
   double y;     /* y component of vector. */
   double z;     /* z component of vector. */
} vector;

void zerovec(vector *r);    /* Set components of vector to zero */
vector vecadd(vector a, vector b);   /* Vector addition a+b. */
vector vecmult(double d, vector a);  /* Multiplication: da. */
vector vecdiv(vector a, double d);  /* Division: a/d. */
vector crossprod(vector a, vector b); /* Cross product, axb */
double magnitude(vector a);  /* Magnitude of vector a */

main(int argc, char *argv[ ])
{

   int i,j,k;
   double j0,j1,j2;
   double m[N_obj],dist1[N_obj],dist2[N_obj];
   double L,minternal,energy;
   double a1[N_obj],a2[N_obj],e1[N_obj],e2[N_obj];
   vector *r,*v,angmom;
   double j10,j11,j12,j13;
   char fp[50];
   FILE *output,*pkd;

   output=fopen("escapecounter.dat","w");
   r=calloc(N_obj,sizeof(vector));
   v=calloc(N_obj,sizeof(vector));
   minternal = M;
   for (k=0; k<N_obj; k++) { 
     dist1[k] = 1.0;
     a1[k] = 1.0;
     e1[k] = 1.0;
   }

   for (i=Out_int; i<=N_steps; i+=Out_int)
    { 
      if (i<10000) sprintf(fp,"ss.0000%d.bt",i);
      else if (i<100000) sprintf(fp,"ss.000%d.bt",i);
      else if (i<1000000) sprintf(fp,"ss.00%d.bt",i);
      else if (i<10000000) sprintf(fp,"ss.0%d.bt",i);
      else sprintf(fp,"ss.%d.bt",i);
      pkd=fopen(fp,"r");
      if (!fp){
	printf("File cannot be opened\n");
	return EXIT_FAILURE;
      }

      for (j=0; j<N_obj; j++)
	{
	  fscanf(pkd,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
		 &j0,&j1,&m[j],&j2,&r[j].x,&r[j].y,&r[j].z,&v[j].x,&v[j].y,&v[j].z,&j10,&j11,&j12,&j13);
	  dist2[j] = magnitude(r[j]);
	}
      fclose(pkd);
      
   for (k=0; k<N_obj; k++) {
     zerovec(&angmom);
     angmom=vecadd(angmom,vecmult(m[k],crossprod(r[k],v[k])));
     L=magnitude(angmom)*2*PI;
     energy=0.5*m[k]*magnitude(v[k])*magnitude(v[k])*4*PI*PI;
     energy-=G*minternal*m[k]/dist2[k];
     a2[k]=-G*minternal*m[k]/(2.0*energy);
     e2[k]=1.0-L*L/(m[k]*m[k]*minternal*G*a2[k]);
     e2[k]=sqrt(e2[k]);
     if (dist1[k]>34000.0 && dist2[k]<20000.0) {
       fprintf(output,"%d %d %lg %lg %lg %lg\n",i,k,dist1[k],dist2[k],a1[k],a1[k]*(1-e1[k]));
     }
     dist1[k] = dist2[k];
     a1[k] = a2[k];
     e1[k] = e2[k];
   }
   
   }
   fclose(output);
   free(r);
   free(v);
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

vector vecmult(double d, vector a)
{
   vector mult;

   mult.x=d*a.x;
   mult.y=d*a.y;
   mult.z=d*a.z;
   return mult;
}

vector vecdiv(vector a, double d)
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

double magnitude(vector a)
{
   double mag=0.0;

   mag=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
   return mag;
}
