/*Set the inital orbital elements for one star.  The rest are far away and circular.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 6.673e-8      /* Newton's constant in cgs */
#define Msun 1.989e33   /* Solar mass in grams */
#define m 1.0      /* Mass of stars in solar masses; assumed constant */
#define alpha 2.0  /* Slope of number density power law */
#define AU 1.496e13  /* AU in cm */
#define amin 100   /* Minimum allowed semimajor axis in AU */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define yr2p (0.5*3.1557e7/PI)  /* One year over 2 pi, in seconds */
#define size1 1.0 /* Size of normal stars for pkdgrav vizualization */
#define size2 1.0 /* Size of binaries */
#define color1 3.0 /* Color of normal stars for pkdgrav vizualization */
#define color2 3.0 /* Color of binaries */

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
  if(argc != 2) {
    printf("Usage: %s Initfile.bt\n", argv[0]);
    return 1;
  }

   char *Initfile = argv[1];
   double M = 1.0e5;
   int N = 1;
   int i,j,iPid;
   double minternal,avar,prob,probmax,xrand;
   double rinfl, amax;
   double E[N],h,R[N],p;
   elements *b;
   vector r,v;
   FILE *initial,*elem;
   
   srand48(iPid=getpid());
   initial=fopen(Initfile,"w");
   elem=fopen("elements.dat","w");

   b=calloc(N,sizeof(elements));

   rinfl = 2.063e5*sqrt(M/1.0e6);
   amax = rinfl*(N/M);
   
   /* create elements array*/
   b[0].a = 500.0;
   b[0].a*=AU;
   b[0].e = 0.99998;
   b[0].nu = 3.14;
   b[0].Om = 2.0*PI*drand48();
   b[0].om = 2.0*PI*drand48();
   b[0].inc = 0.807048;
   E[0]=acos((b[0].e+cos(b[0].nu))/(1.0+b[0].e*cos(b[0].nu)));
   if (b[0].nu>PI) E[0]=2.0*PI-E[0];
   R[0]=b[0].a*(1.0-b[0].e*cos(E[0]));
   fprintf(elem,"%d %lg %lg %lg %lg %lg %lg %lg\n",0,b[0].a/AU,b[0].e,b[0].inc,b[0].Om,b[0].om,b[0].nu,R[0]/AU);

   for (i=1; i<N; i++)
   {
      avar=drand48();
      b[i].a=3.0e4+100.0*avar;
      b[i].e=0;
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

      fprintf(elem,"%d %lg %lg %lg %lg %lg %lg\n",i,b[i].a/AU,b[i].e,b[i].inc,b[i].Om,b[i].om,b[i].nu);

   }

   for (i=0; i<N; i++)
     {
       /*calculate minternal*/
       minternal = M*Msun;
       for (j=0; j<N; j++)
	 {
	   if (R[i] > R[j]) {
	     minternal+=m*Msun; //add mass of stars
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
       /*convert to pkdgrav units*/
       r.x/=AU;
       r.y/=AU;
       r.z/=AU;
       v.x*=yr2p/AU;
       v.y*=yr2p/AU;
       v.z*=yr2p/AU;

       fprintf(initial,"%i %i %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
	       i,i,m,size1,r.x,r.y,r.z,v.x,v.y,v.z,0.0,0.0,0.0,color1);

     }


   fclose(initial);
   fclose(elem);
}
