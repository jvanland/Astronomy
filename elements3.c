/* Prints out the orbital elements for a
collection of objects. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.673e-8      /* Newton's constant in cgs */
#define Msun 1.989e33   /* Solar mass in grams */
#define M 1.0e5    /* Mass of SMBH in solar masses */
#define N 10000     /* Number of lines */
#define AU 1.496e13  /* AU in cm */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define yr2p (0.5*3.1557e7/PI)  /* One year over 2 pi, in seconds */


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
   int i,t;
   double L,m0,r2,v2,lx,ly,lz,sma,l2,ecc,peri;
   float h1,h2,h4,h11,h12,h13,h14;
   vector r,v;
   double energy;
   double a,e,incl;
   vector angmom,amomhat;
   FILE *initial;

   initial=fopen(Initfile,"r");
   zerovec(&amomhat);
   zerovec(&angmom);
   for (i=0; i<N; i++)
   {
     fscanf(initial,"%f %f %lg %f %lg %lg %lg %lg %lg %lg %f %f %f %f",
		 &h1,&h2,&m0,&h4,&r.x,&r.y,&r.z,&v.x,&v.y,&v.z,&h11,&h12,&h13,&h14);
      angmom=vecmult(m0,crossprod(r,v));
      amomhat = vecdiv(angmom,magnitude(angmom));
      incl=acos(amomhat.z);
      L=magnitude(angmom)*Msun*AU*AU/yr2p;
      energy=0.5*m0*Msun*magnitude(v)*magnitude(v)*AU*AU/(yr2p*yr2p);
      energy-=G*M*Msun*m0*Msun/(magnitude(r)*AU);
      a=-G*M*Msun*m0*Msun/(2.0*energy)/AU;
      e=1.0-L*L/(Msun*m0*m0*Msun*Msun*M*G*a*AU);
      e=sqrt(e);
      r2 = r.x*r.x + r.y*r.y + r.z*r.z;
      v2 = v.x*v.x + v.y*v.y + v.z*v.z;
      energy = 0.5*v2 - M/sqrt(r2);
      sma = -0.5*M/energy;
      lx = r.y*v.z - r.z*v.y;
      ly = r.z*v.x - r.x*v.z;
      lz = r.x*v.y - r.y*v.x;
      l2 = lx*lx+ly*ly+lz*lz;
      ecc = sqrt(1.0-l2/(M*sma));
      peri = sma*(1.0-ecc);
      printf("%d %lg %lg %lg %lg %lg %lg\n",i,a,e,a*(1.0-e),sma,ecc,peri);
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

