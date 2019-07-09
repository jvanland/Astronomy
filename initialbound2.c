/* Draw a binary of 10 Msun BH's with a = 1 AU, the 
eccentricity from a thermal distribution P(e)=2e, 
the true anomaly using a proper time weighting,
and the other angles from an isotropic distribution.
These are then changed to Cartesian coordinates for the
output. 
Called by driver.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 39.47841760435743 /* Newton's constant for AU,yr,Msun*/
#define M 10.0    /* Mass of the more massive member of the binary */
#define m 10.0      /* Mass of the less massive member of the binary */
#define N 1      /* Number of binaries */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define yr2p (0.5/PI)  /* One year over 2 pi, in years */

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

void convert(double m1, double m2, elements b, vector *r, vector *v);

main()
{
   int i,iPid;
   double m1,m2,prob,probmax,xrand;
   elements b;
   vector r,v;
   FILE *initial;
   
   srand48(iPid=getpid());
   initial=fopen("initial0.dat","w");

   for (i=0; i<N; i++)
   {
      b.a=1.0;  /* Fix binary semimajor axis at 1 AU */
      m1=M;
      m2=m;
      b.e=sqrt(drand48());
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
      convert(m1,m2,b,&r,&v);
      v.x*=yr2p;
      v.y*=yr2p;
      v.z*=yr2p;
      r.x*=m2/(m1+m2);
      r.y*=m2/(m1+m2);
      r.z*=m2/(m1+m2);
      v.x*=m2/(m1+m2);
      v.y*=m2/(m1+m2);
      v.z*=m2/(m1+m2);
      fprintf(initial,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
         m1,r.x,r.y,r.z,v.x,v.y,v.z);
      fprintf(initial,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
         m2,-(m1/m2)*r.x,-(m1/m2)*r.y,-(m1/m2)*r.z,
         -(m1/m2)*v.x,-(m1/m2)*v.y,-(m1/m2)*v.z);
   }
   fclose(initial);
}

void convert(double m1, double m2, elements b, vector *r, vector *v)
{
/* Converts from orbital elements to vector position and velocity. */

   double h,R,p,Mtot;

   Mtot=m1+m2;

/* First the eccentric anomaly. */

   R=b.a*(1.0-b.e*b.e)/(1+b.e*cos(b.nu));
   h=sqrt(G*Mtot*b.a*(1.0-b.e*b.e));
   p=b.a*(1.0-b.e*b.e);
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
