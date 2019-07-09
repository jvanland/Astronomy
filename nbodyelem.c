/* Converts Cartesian coordinates to Orbital elements 
for a selection of objects in bt format. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 1.0 /* Newton's constant for AU,yr/2pi,Msun */
#define M 1.0e3    /* Mass of SMBH in solar masses */
#define N 1001     /* Number of stars */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define Rconv 0.000440724
#define Vconv 1.0601
#define Tconv 0.000415736
#define Mtot 2019


typedef struct
{
   double x;     /* x component of vector. */
   double y;     /* y component of vector. */
   double z;     /* z component of vector. */
} vector;

void zerovec(vector *r);    /* Set components of vector to zero */
vector vecadd(vector a, vector b);   /* Vector addition a+b. */
vector vecsub(vector a, vector b);   /* Vector subtraction a-b. */
vector vecmult(double d, vector a);  /* Multiplication: da. */
vector vecdiv(vector a, double d);  /* Division: a/d. */
vector crossprod(vector a, vector b); /* Cross product, axb */
double magnitude(vector a);  /* Magnitude of vector a */

int main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Initfile\n", argv[0]);
    return 1;
  }
   char *Initfile = argv[1];
   int i,j;
   double m[N],L,dist[N], minternal;
   vector *r,*v,*rcent,*vcent;
   double a,e,incl,energy;
   vector angmom,amomhat;
   FILE *initial;

   initial=fopen(Initfile,"r");
   r=calloc(N,sizeof(vector));
   v=calloc(N,sizeof(vector));
   rcent=calloc(N,sizeof(vector));
   vcent=calloc(N,sizeof(vector));
   zerovec(&angmom);
   for (i=0; i<N; i++)
   {
      fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg",
	     &m[i],&r[i].x,&r[i].y,&r[i].z,&v[i].x,&v[i].y,&v[i].z);
      rcent[i] = vecsub(r[i],r[0]);
      vcent[i] = vecsub(v[i],v[0]);
      dist[i] = magnitude(rcent[i]);
   }
   
   for (i=0; i<N; i++)
     {
       minternal = 0.0;
       for (j=0; j<N; j++)
	 {
	   if (dist[i] > dist[j]) {
	     minternal+=m[j];
	   }
	 }
      zerovec(&angmom);
      angmom=vecadd(angmom,vecmult(m[i],crossprod(rcent[i],vcent[i])));
      amomhat = vecdiv(angmom,magnitude(angmom));
      incl=acos(amomhat.z)*180.0/PI;
      L=magnitude(angmom);
      energy=0.5*m[i]*magnitude(vcent[i])*magnitude(vcent[i]);
      energy-=G*minternal*m[i]/dist[i];
      a=-G*minternal*m[i]/(2.0*energy);
      e=1.0-L*L/(m[i]*m[i]*minternal*G*a);
      e=sqrt(e);
      /* Conversions */
      dist[i] = dist[i]/Rconv;
      a = a/Rconv;
      printf("%d %lg %lg %lg %lg %lg\n",i,a,e,incl,minternal,dist[i]);
   }
   fclose(initial);
   free(r);
   free(v);
   free(rcent);
   free(vcent);
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
