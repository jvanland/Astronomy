/* Prints out the orbital elements for a
collection of objects. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 39.47841760435743 /* Newton's constant for AU,yr,Msun */
#define M 5.0e3    /* Mass of SMBH in solar masses */
#define m 20.0    /* mass of binary particle */
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

int main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Initfile\n", argv[0]);
    return 1;
  }
   char *Initfile = argv[1];
   int i,t,N=0;
   double L,m0;
   vector r,v;
   double energy,dist;
   double a,e,incl;
   vector angmom,amomhat;
   FILE *initial;
   char buf[5000];

   initial=fopen(Initfile,"r");
   while((fgets(buf,5000,initial)) != NULL) {
     N++;
   }
   rewind(initial);
   zerovec(&amomhat);
   zerovec(&angmom);
   for (i=0; i<N; i++)
   {
     fscanf(initial,"%d %lg %lg %lg %lg %lg %lg",&t,&r.x,&r.y,&r.z,&v.x,&v.y,&v.z);
     dist = magnitude(r);
     angmom=vecmult(m,crossprod(r,v));
     amomhat = vecdiv(angmom,magnitude(angmom));
     incl=acos(amomhat.z)*180.0/PI;
     L=magnitude(angmom)*2*PI;
     energy=0.5*m*magnitude(v)*magnitude(v)*4*PI*PI;
     energy-=G*M*m/dist;
     a=-G*M*m/(2.0*energy);
     e=1.0-L*L/(m*m*M*G*a);
     e=sqrt(e);
     printf("%d %lg %lg %lg %lg %lg\n",i,a,e,incl,a*(1-e),dist);
   }
   fclose(initial);
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

