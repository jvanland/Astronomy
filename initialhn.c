/* New approach to initial conditions file.  We generate
orbital elements for the stars around an SMBH, drawing
the semimajor axis from a power law distribution with
index -alpha, the eccentricity from a thermal distribution
P(e)=2e, the true anomaly using a proper time weighting,
and the other angles from an isotropic distribution.
These are then changed to Cartesian coordinates for the
output. Output is formatted for hnbody.  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 6.673e-8      /* Newton's constant in cgs */
#define Msun 1.989e33   /* Solar mass in grams */
#define m 1.0      /* Mass of stars in solar masses; assumed constant */
#define mbin 20.0    /*Mass of BH binary and mirror */
#define alpha 2.0  /* Slope of number density power law */
#define AU 1.496e13  /* AU in cm */
#define amin 100   /* Minimum allowed semimajor axis in AU */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define yr 3.1557e7  /* One year, in seconds */

typedef struct
{
   double x;     /* x component of vector. */
   double y;     /* y component of vector. */
   double z;     /* z component of vector. */
} vector;

typedef struct
{
   double a;     /* Semimajor axis. */
   double e;     /* Eccentricity. */
   double inc;     /* Inclination. */
   double Om;    /* Ascending node. */
   double om;    /* Argument of pericenter. */
   double nu;    /* True anomaly. */
} elements;

main(int argc, char *argv[ ])
{
  if(argc != 3) {
    printf("Usage: %s MSMBH(Solar Masses, dbl) NBin(int)\n", argv[0]);
    return 1;
  }

   double M = atof(argv[1]);
   int N = M/10;
   int NBin = atoi(argv[2]);
   int i,j,iPid;
   double minternal,avar,prob,probmax,xrand;
   double rinfl, amax;
   double E[N],h,R[N],p;
   elements *b;
   vector r,v;
   FILE *binary,*star;
   
   srand48(10);
   binary=fopen("binaries.hnb","w");
   star=fopen("stars.hnb","w");

   b=calloc(N,sizeof(elements));

   rinfl = 2.063e5*sqrt(M/1.0e6);
   amax = rinfl*(N/M);
   
   /* create elements array*/
   for (i=0; i<N; i++)
   {
      avar=drand48();
      do {
	b[i].a=pow((pow(amin,3.0-alpha)+avar*(pow(amax,3.0-alpha)-pow(amin,3.0-alpha))),1.0/(3.0-alpha));
	b[i].e=sqrt(drand48());
      } while (b[i].a*(1.0-b[i].e)<10.0);
      b[i].a*=AU;
      probmax=(1.0+b[i].e)*(1.0+b[i].e);
      do
      {
         b[i].nu=2.0*PI*drand48();
         prob=(1.0-b[i].e*b[i].e)/(1.0+b[i].e*cos(b[i].nu));
         prob=prob*prob;
         xrand=drand48();
      } while (xrand>prob/probmax);
      b[i].Om=2.0*PI*drand48();
      b[i].om=2.0*PI*drand48();
      b[i].inc=acos(1.0-2.0*drand48());
      if (drand48()<0.5) b[i].inc*=-1.0;

      /*calculate eccentric anomaly and distance from SMBH for each star*/
      E[i]=acos((b[i].e+cos(b[i].nu))/(1.0+b[i].e*cos(b[i].nu)));
      if (b[i].nu>PI) E[i]=2.0*PI-E[i];
      R[i]=b[i].a*(1.0-b[i].e*cos(E[i]));

   }

   for (i=0; i<N; i++)
     {
       /*calculate minternal*/
       minternal = M*Msun;
       if (i >= N-NBin) { // stars are LWP and only feel SMBH
	 for (j=N-NBin; j<N; j++)
	   {
	     if (R[i] > R[j]) {
	       minternal+=mbin*Msun; //add mass of binaries if interior
	     }       
	   }
	 for (j=0; j<N-NBin; j++)
	   {
	     if (R[i] > R[j]) {
	       minternal+=m*Msun; //add mass of stars
	     }
	   }
       }
       /*convert to cartesian coordinates*/
       h=sqrt(G*minternal*b[i].a*(1.0-b[i].e*b[i].e));
       p=(1.0+b[i].e*cos(b[i].nu))*b[i].a*(1.0-b[i].e*cos(E[i]));
       
       r.x=R[i]*(cos(b[i].Om)*cos(b[i].om+b[i].nu)-sin(b[i].Om)*sin(b[i].om+b[i].nu)*cos(b[i].inc));
       r.y=R[i]*(sin(b[i].Om)*cos(b[i].om+b[i].nu)+cos(b[i].Om)*sin(b[i].om+b[i].nu)*cos(b[i].inc));
       r.z=R[i]*(sin(b[i].inc)*sin(b[i].om+b[i].nu));

       v.x= r.x*h*b[i].e*sin(b[i].nu)/(R[i]*p);
       v.x-=(h/R[i])*(cos(b[i].Om)*sin(b[i].om+b[i].nu)+sin(b[i].Om)*cos(b[i].om+b[i].nu)*cos(b[i].inc));
       v.y= r.y*h*b[i].e*sin(b[i].nu)/(R[i]*p);
       v.y-=(h/R[i])*(sin(b[i].Om)*sin(b[i].om+b[i].nu)-cos(b[i].Om)*cos(b[i].om+b[i].nu)*cos(b[i].inc));
       v.z= r.z*h*b[i].e*sin(b[i].nu)/(R[i]*p);
       v.z+=(h/R[i])*sin(b[i].inc)*cos(b[i].om+b[i].nu);
       /*convert to hnbody units*/
       r.x/=AU;
       r.y/=AU;
       r.z/=AU;
       v.x*=yr/AU;
       v.y*=yr/AU;
       v.z*=yr/AU;

       if (i < N-NBin) {
	 fprintf(star,"%lg %lg %lg %lg %lg %lg %lg\n",m,r.x,r.y,r.z,v.x,v.y,v.z);
       }
       /* The last stars are the binaries */
       else {
	 fprintf(binary,"%lg %lg %lg %lg %lg %lg %lg\n",mbin,r.x,r.y,r.z,v.x,v.y,v.z);
       }

     }


   fclose(star);
   fclose(binary);
}
