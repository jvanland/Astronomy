/* Converts Cartesian coordinates to Orbital elements 
for a selection of objects in bt format. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 39.47841760435743 /* Newton's constant for AU,yr,Msun */
#define M 1.0e3    /* Mass of SMBH in solar masses */
#define N 1000     /* Number of stars */
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
  if(argc != 2) {
    printf("Usage: %s Initfile\n", argv[0]);
    return 1;
  }
   char *Initfile = argv[1];
   int i,j;
   double m[N],L,dist[N], minternal;
   vector *r,*v;
   double a,e,incl,energy;
   vector angmom,amomhat;
   FILE *initial;

   initial=fopen(Initfile,"r");
   r=calloc(N,sizeof(vector));
   v=calloc(N,sizeof(vector));
   zerovec(&angmom);
   for (i=0; i<N; i++)
   {
      fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg",
	     &m[i],&r[i].x,&r[i].y,&r[i].z,&v[i].x,&v[i].y,&v[i].z);
      dist[i] = magnitude(r[i]);
   }
   
   for (i=0; i<N; i++)
     {
       minternal = M;
       for (j=0; j<N; j++)
	 {
	   if (dist[i] > dist[j]) {
	     minternal+=m[j];
	   }
	 }
      zerovec(&angmom);
      angmom=vecadd(angmom,vecmult(m[i],crossprod(r[i],v[i])));
      amomhat = vecdiv(angmom,magnitude(angmom));
      incl=acos(amomhat.z)*180.0/PI;
      L=magnitude(angmom)*2*PI;
      energy=0.5*m[i]*magnitude(v[i])*magnitude(v[i])*4*PI*PI;
      energy-=G*minternal*m[i]/dist[i];
      a=-G*minternal*m[i]/(2.0*energy);
      e=1.0-L*L/(m[i]*m[i]*minternal*G*a);
      e=sqrt(e);
      printf("%d %lg %lg %lg %lg %lg\n",i,a,e,incl,a*(1-e),dist[i]);
   }
   fclose(initial);
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
