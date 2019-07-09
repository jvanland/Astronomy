/* New approach to initial conditions file.  We generate
orbital elements for the stars around an SMBH, drawing
the semimajor axis from a power law distribution with
index -alpha, the eccentricity from a thermal distribution
P(e)=2e, the true anomaly using a proper time weighting,
and the other angles from an isotropic distribution.
These are then changed to Cartesian coordinates for the
output. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 6.673e-8      /* Newton's constant in cgs */
#define Msun 1.989e33   /* Solar mass in grams */
#define M 1.0e4    /* Mass of SMBH in solar masses */
#define N 1000     /* Number of stars */
#define m 1.0      /* Mass of stars in solar masses; assumed constant */
#define alpha 2.0  /* Slope of number density power law */
#define AU 1.496e13  /* AU in cm */
#define amin 100   /* Minimum allowed semimajor axis in AU */
#define amax 2.1e3   /* Maximum allowed semimajor axis in AU */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define yr2p (0.5*3.1557e7/PI)  /* One year over 2 pi, in seconds */

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
   double i;     /* Inclination. */
   double Om;    /* Ascending node. */
   double om;    /* Argument of pericenter. */
   double nu;    /* True anomaly. */
} elements;

void convert(double Mtot, elements b, vector *r, vector *v);

main()
{
   int i,iPid;
   double mtot,avar,prob,probmax,xrand;
   elements b;
   vector r,v;
   FILE *initial, *elements;
   
   srand48(iPid=getpid());
   initial=fopen("init.dat","w");
   elements=fopen("elements.dat","w");

   for (i=0; i<N/2; i++)
   {
      avar=drand48();
      do {
	b.a=pow((pow(amin,3.0-alpha)+avar*(pow(amax,3.0-alpha)-pow(amin,3.0-alpha))),1.0/(3.0-alpha));
	b.e=sqrt(drand48());
	//b.e=0;
      } while (b.a*(1.0-b.e)<10.0);
      b.a*=AU;
      mtot=M+m*N*avar;
      mtot*=Msun;
      probmax=(1.0+b.e)*(1.0+b.e);
      do
      {
         b.nu=2.0*PI*drand48();
         prob=(1.0-b.e*b.e)/(1.0+b.e*cos(b.nu));
         prob=prob*prob;
         xrand=drand48();
      } while (xrand>prob/probmax);
      b.Om=2.0*PI*drand48();
      b.om=2.0*PI*drand48();
      b.i=acos(1.0-2.0*drand48());
      if (drand48()<0.5) b.i*=-1.0;
      convert(mtot,b,&r,&v);
      r.x/=AU;
      r.y/=AU;
      r.z/=AU;
      v.x*=yr2p/AU;
      v.y*=yr2p/AU;
      v.z*=yr2p/AU;
      fprintf(initial,"%lg %lg %lg %lg %lg %lg %lg\n",
         m,r.x,r.y,r.z,v.x,v.y,v.z);
      fprintf(initial,"%lg %lg %lg %lg %lg %lg %lg\n",
         m,-r.x,-r.y,-r.z,-v.x,-v.y,-v.z);
      fprintf(elements,"%d %lg %lg %lg\n",i,b.a/AU,b.e,b.i);
   }
   fclose(initial);
   fclose(elements);
}

void convert(double Mtot, elements b, vector *r, vector *v)
{
/* Converts from orbital elements to vector position and velocity. */

   double E,h,R,p;

/* First the eccentric anomaly. */

   E=acos((b.e+cos(b.nu))/(1.0+b.e*cos(b.nu)));
   if (b.nu>PI) E=2.0*PI-E;

   R=b.a*(1.0-b.e*cos(E));
   h=sqrt(G*Mtot*b.a*(1.0-b.e*b.e));
   p=(1.0+b.e*cos(b.nu))*b.a*(1.0-b.e*cos(E));
   (*r).x=R*(cos(b.Om)*cos(b.om+b.nu)-sin(b.Om)*sin(b.om+b.nu)*cos(b.i));
   (*r).y=R*(sin(b.Om)*cos(b.om+b.nu)+cos(b.Om)*sin(b.om+b.nu)*cos(b.i));
   (*r).z=R*(sin(b.i)*sin(b.om+b.nu));

   (*v).x= (*r).x*h*b.e*sin(b.nu)/(R*p);
   (*v).x-=(h/R)*(cos(b.Om)*sin(b.om+b.nu)+sin(b.Om)*cos(b.om+b.nu)*cos(b.i));
   (*v).y= (*r).y*h*b.e*sin(b.nu)/(R*p);
   (*v).y-=(h/R)*(sin(b.Om)*sin(b.om+b.nu)-cos(b.Om)*cos(b.om+b.nu)*cos(b.i));
   (*v).z= (*r).z*h*b.e*sin(b.nu)/(R*p);
   (*v).z+=(h/R)*sin(b.i)*cos(b.om+b.nu);
}
