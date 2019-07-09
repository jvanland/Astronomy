/* New approach to initial conditions file.  We generate
orbital elements for the stars around an SMBH, drawing
the semimajor axis from a power law distribution with
index -alpha, the eccentricity from a thermal distribution
P(e)=2e, the true anomaly using a proper time weighting,
and the other angles from an isotropic distribution.
These are then changed to Cartesian coordinates for the
output. To force the total angular momentum to be zero, 
each star has a mirror.  A BH binary particle (and mirror) 
is also included.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 1.0 /* Newton's constant for AU,yr/2pi,Msun */
#define m 1.0      /* Mass of stars in solar masses; assumed constant */
#define mbin 20.0    /* Mass of binary in solar masses */
#define alpha 2.0  /* Slope of number density power law */
#define amin 10   /* Minimum allowed semimajor axis in AU */
#define PI 3.141592653589793  /* pi to unnecessary precision */

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

int main(int argc, char *argv[ ])
{
  if(argc != 4) {
    printf("Usage: %s MSMBH(Solar Masses, dbl) Nstars NBin(int)\n", argv[0]);
    return 1;
  }

   double M = atof(argv[1]);
   int N = atoi(argv[2]);
   int NBin = atoi(argv[3]);
   int i,j,cnt1,cnt2,iPid;
   double minternal,prob,probmax,xrand;
   double rinfl,amax,Lrel;
   double E[N],h,R[N],p;
   elements *b;
   vector r,v;
   FILE *initial,*elem;
   
   srand48(iPid=getpid());
   initial=fopen("init.dat","w");
   elem=fopen("elements.dat","w");
   fprintf(initial,"%lg 0.0 0.0 0.0 0.0 0.0 0.0\n", M);

   b=calloc(N,sizeof(elements));

   rinfl = 2.063e5*sqrt(M/1.0e6);
   amax = rinfl*(N/M);
   //  amax = rinfl*0.307;   

   /* create elements array*/
   for (i=0; i<N; i++)
   {
     //      do {
	b[i].a=pow(drand48()*(pow(amax,3.0-alpha)),1.0/(3.0-alpha));
	b[i].e=sqrt(drand48());
	if(i >= N-NBin) {
	  //	  do {
	  b[i].a=1500.0;
	    //	    b[i].e = sqrt(drand48());
	    //	  } while (b[i].a*(1-b[i].e) < 0.5*b[i].a);
	}
	Lrel = (1.0-b[i].e*b[i].e)*sqrt(b[i].a/amax)*(mbin/m);
	//      } while (Lrel < 1.0);
	if(Lrel < 1.0) cnt1 ++;
	if(Lrel > 1.0) cnt2 ++;
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

      fprintf(elem,"%d %lg %lg %lg %lg %lg\n",i,b[i].a,b[i].e,b[i].inc,Lrel,R[i]);

   }
   printf("Lrel < 1.0 = %d Lrel > 1.0 = %d\n",cnt1,cnt2);

   for (i=0; i<N; i++)
     {
       /*calculate minternal*/
       minternal = M;
       for (j=N-NBin-1; j<N; j++) 
         {
	   if (R[i] > R[j]) {
	     minternal+=mbin; //add mass of binaries if interior
           }       
         }
       for (j=0; j<N-NBin; j++)
	 {
	   if (R[i] > R[j]) {
	     minternal+=m; //add mass of stars
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

       if (i < N-NBin) {
	 fprintf(initial,"%lg %lg %lg %lg %lg %lg %lg\n",
		 m,r.x,r.y,r.z,v.x,v.y,v.z);
       }

       /* The last stars are the binaries */
       else {
	 fprintf(initial, "%lg %lg %lg %lg %lg %lg %lg\n",
		 mbin,r.x,r.y,r.z,v.x,v.y,v.z);
       }
     }

   fclose(initial);
   fclose(elem);
}
