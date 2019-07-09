/* 
Draw a binary around a central SMBH.
For use in hndrag.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 39.47841760435743 /* Newton's constant for AU,yr,Msun*/
#define M 10.0    /* Mass of the more massive member of the binary */
#define m 10.0      /* Mass of the less massive member of the binary */
#define M_SMBH 5.0e3   /* Mass of the SMBH */
#define N 1      /* Number of binaries */
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
   double i;     /* Inclination. */
   double Om;    /* Ascending node. */
   double om;    /* Argument of pericenter. */
   double nu;    /* True anomaly. */
} elements;

void convert(double m1, double m2, elements b, vector *r, vector *v);

int main()
{
   int i,iPid;
   double m1,m2,mbin,msmbh;
   elements b1,b2;
   vector r1,r2,rsuper,v1,v2,vsuper;
   FILE *initial;
   
   initial=fopen("binary.hnb","w");
   /* Define orbit of binary */
   b1.a=0.01;
   b1.e=0.0;
   b1.i=90.0*PI/180.0;
   b1.Om=0.0;
   b1.om=0.0;
   b1.nu=0.0;
   m1=M;
   m2=m;
   /* Convert to cartesian coordinates and put c.o.m. at 0 */
   convert(m1,m2,b1,&r1,&v1);
   r1.x*=m2/(m1+m2);
   r2.x=-(m1/m2)*r1.x;
   r1.y*=m2/(m1+m2);
   r2.y=-(m1/m2)*r1.y;
   r1.z*=m2/(m1+m2);
   r2.z=-(m1/m2)*r1.z;
   v1.x*=m2/(m1+m2);
   v2.x=-(m1/m2)*v1.x;
   v1.y*=m2/(m1+m2);
   v2.y=-(m1/m2)*v1.y;
   v1.z*=m2/(m1+m2);
   v2.z=-(m1/m2)*v1.z;
   /* Define super-orbit */
   b2.a=1500.0;
   b2.e=0.0;
   b2.i=0.0;
   b2.Om=0.0;
   b2.om=0.0;
   b2.nu=0.0;
   mbin=m1+m2;
   msmbh=M_SMBH;
   /* Convert to cartesian coordinates */
   convert(mbin,msmbh,b2,&rsuper,&vsuper);
   /* Place binary orbit on top of super-orbit */
   r1.x+=rsuper.x;
   r2.x+=rsuper.x;
   r1.y+=rsuper.y;
   r2.y+=rsuper.y;
   r1.z+=rsuper.z;
   r2.z+=rsuper.z;
   v1.x+=vsuper.x;
   v2.x+=vsuper.x;
   v1.y+=vsuper.y;
   v2.y+=vsuper.y;
   v1.z+=vsuper.z;
   v2.z+=vsuper.z;

   fprintf(initial,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
	   m1,r1.x,r1.y,r1.z,v1.x,v1.y,v1.z);
   fprintf(initial,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
	   m2,r2.x,r2.y,r2.z,v2.x,v2.y,v2.z);

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
