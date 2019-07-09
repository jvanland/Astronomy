/* Prints out the net momentum and angular momentum for a
collection of objects. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.673e-8      /* Newton's constant in cgs */
#define Msun 1.989e33   /* Solar mass in grams */
#define M 1.0e4    /* Mass of SMBH in solar masses */
#define N 1000      /* Number of stars */
#define m0 1.0      /* Mass of stars in solar masses; assumed constant */
#define AU 1.496e13  /* AU in cm */
#define amax 2.1e3 /* Maximum allowed semimajor axis in AU */
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
vector crossprod(vector a, vector b); /* Cross product, axb */
double magnitude(vector a);  /* Magnitude of vector a */

main()
{
   int i;
   double j0,j1,j2;
   double m[N],L;
   vector *r,*v;
   double j10,j11,j12,j13;
   vector mom,angmom;
   FILE *initial;

   initial=fopen("ss.90000.bt","r");
   r=calloc(N,sizeof(vector));
   v=calloc(N,sizeof(vector));
   zerovec(&mom);
   zerovec(&angmom);
   for (i=0; i<N; i++)
   {
      fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
	     &j0,&j1,&m[i],&j2,&r[i].x,&r[i].y,&r[i].z,&v[i].x,&v[i].y,&v[i].z,&j10,&j11,&j12,&j13);
      mom=vecadd(mom,vecmult(m[i],v[i]));
      angmom=vecadd(angmom,vecmult(m[i],crossprod(r[i],v[i])));
   }
   fclose(initial);
   printf("Total net momentum is %lg x + %lg y + %lg z\n",
      mom.x,mom.y,mom.z);
   printf("Total net ang mom is %lg x + %lg y + %lg z\n",
      angmom.x,angmom.y,angmom.z);
   L=magnitude(angmom);
   printf("Angular momentum is equivalent to %lg stars of mass %lg Msun\n",
      L*Msun*AU*AU/yr2p/(m0*Msun*sqrt(G*M*Msun*amax*AU)),m0);
   printf("in a circular orbit at %lg AU.\n",amax);
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
