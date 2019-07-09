/* Initial condition generator.  We assume that a supermassive
black hole of mass M solar masses has N stars around it of
mass m.  These stars have a number density n~r^{-alpha}, e.g.,
alpha=7/4 gives the Bahcall-Wolf distribution.  We assume that
the stars are distributed in a spherically symmetric way, and
that the velocity distribution is isotropic at all points, with
a magnitude drawn from a Maxwell-Boltzmann distribution with
a mean equal to the velocity of a circular orbit at that
location, including the SMBH mass and the masses of the stars
interior to the orbit (this guarantees that the virial theorem
holds).  Output is in units of AU for distance, and AU/yr for
the three velocity components. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 6.673e-8      /* Newton's constant in cgs */
#define Msun 1.989e33   /* Solar mass in grams */
#define M 1.0e6    /* Mass of SMBH in solar masses */
#define N 10000      /* Number of stars */
#define m 1.0      /* Mass of stars in solar masses; assumed constant */
#define alpha 7/4  /* Slope of number density power law */
#define AU 1.496e13  /* AU in cm */
#define rmin 100   /* Minimum allowed radius in AU */
#define rmax 2.0e5 /* Maximum allowed radius in AU */
#define yr 3.1557e7  /* One year, in seconds */
#define PI 3.141592653589793  /* pi to unnecessary precision */

double gasdev(void);

main()
{
   int i,iPid;
   double rvar;
   double r,ctheta,stheta,phi,x,y,z,mint,v1d,vx,vy,vz;
   FILE *initial;
   
   srand48(iPid=getpid());
   initial=fopen("initial.dat","w");

   for (i=0; i<N; i++)
   {
      ctheta=1.0-2.0*drand48();
      stheta=sqrt(1.0-ctheta*ctheta);
      phi=2.0*PI*drand48();
      rvar=drand48();
      r=pow((pow(rmin,3.0-alpha)+rvar*(pow(rmax,3.0-alpha)-pow(rmin,3.0-alpha))),1.0/(3.0-alpha));
      x=r*stheta*cos(phi);
      y=r*stheta*sin(phi);
      z=r*ctheta;
/* x, y, and z are in AU now */
      mint=M+m*N*rvar;
      v1d=(1.0/sqrt(3.0))*sqrt(G*mint*Msun/(r*AU));  /* 1-D vel disp in cm/s */
      v1d*=yr/AU;   /* And in AU/yr */
      vx=v1d*gasdev();
      vy=v1d*gasdev();
      vz=v1d*gasdev();
      fprintf(initial,"%lg %lg %lg %lg %lg %lg %lg\n",
         m,x,y,z,vx,vy,vz);
   }
   fclose(initial);
}

double gasdev(void)
{
/* Generate gaussian random deviate with mean zero and variance 1. */

   static double gset;
   double fac,r,v1,v2;

     do
     {
       v1=2.0*drand48()-1.0;
       v2=2.0*drand48()-1.0;
       r=v1*v1+v2*v2;
     } while (r>=1.0 || r==0.0);
     fac=sqrt(-2.0*log(r)/r);
     gset=v1*fac;
     return v2*fac;
}
