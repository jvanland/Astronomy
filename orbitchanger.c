/* Sets up an orbit around a SMBH where
one or more of the orbital elements changes
predicatbly over time.  Prints out the
cartesian coordinates of the orbit.  For
use with hnbody to examine Kozai. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "delaunay.h"
#include "helio.h"

#define G 39.47841760435743  /* Gravitational constant for Msun, AU, yr */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define Mstar 1.0  /* Mass of star */
#define Mbin 20.0  /* Mass of binary */
#define Mbh 1.0e5  /* Mass of SMBH */
#define delta 100.0  /* Time between outputs (yr) */
#define Tfin 1.0e5  /* Total time */


void deltohel(double,struct delaunay *,struct helio *);

main()
{
  int i,N,iPid;
   struct helio hel;
   struct delaunay del;
   double deltaM,deltaA,deltaA2,deltaE,deltaI,deltaOm,deltaw;
   double rinfl,clmb,Nstar;
   double tRR,tRRv,torb,tprec,t=0;
   FILE *orbit,*init,*elem;

   srand48(iPid=getpid());
   rinfl = 2.06265e5*sqrt(Mbh/1.0e6);  /* radius of influence of the SMBH */
   clmb = 10.0;  /* Cuolomb logarithm */

   orbit=fopen("orbit.dat","w");
   elem=fopen("delinit.dat","w");

   /*Read in initial orbital elements */
   init=fopen("eleminit.dat","r");
   fscanf(init,"%lg %lg %lg %lg %lg %lg",
      &del.sma,&del.ecc,&del.inc,&del.lan,&del.lop,&del.mea);
   fclose(init);

   /* Convert to cartesian coordinates */
   deltohel(G*(Mbh+Mbin),&del,&hel);
   fprintf(orbit,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
	   t,hel.x,hel.y,hel.z,hel.vx,hel.vy,hel.vz);

   N = Tfin/delta;

   for (i=0; i<N; i++)
     {
       t += delta;  /* update time */
       deltaA = sqrt(G*del.sma/Mbh)*Mbin*Mstar*clmb*delta/rinfl;
       del.sma -= deltaA;  /* sma shrinks due to dynamical friction */
       //Nstar = Mbh*del.sma/rinfl;  /* Number of stars interior */
       torb = 2*PI*sqrt(del.sma*del.sma*del.sma/(G*Mbh));  /* Orbital period */
       //tRR = (Mbh/Mbin)*torb;  /* Resonant relaxation time */
       //tRRv = sqrt(Nstar)*((Nstar*Mstar+Mbh)/(Nstar*Mstar))*torb;  /* Vector resonant relaxation time */
       //tprec = (Mbh/(Nstar*Mstar))*torb;  /* Precession time */
       //deltaA2 = del.sma*sqrt(delta/tRR);  /* sma also random walks on tRR */
       //if (del.sma < 100.0) del.sma +=deltaA2;
       //else del.sma += (2.0*drand48() - 1.0)*deltaA2;
       //deltaE = (1-del.ecc)*sqrt(delta/tRR)/3.0;
       //if (del.ecc-deltaE < 0.0) del.ecc += deltaE;  /* eccentricity can't go below zero */ 
       //else if (del.ecc+deltaE > 0.999) del.ecc -= deltaE;  /* or above 1 */ 
       //else del.ecc += (2.0*drand48() - 1.0)*deltaE;  /* ecc random walks on tRR */
       //deltaI = PI*sqrt(delta/tRRv)/3.0;
       //del.inc += (2.0*drand48() - 1.0)*deltaI;  /* inc random walks on tRRv */
       //del.inc = fmod(del.inc,2*PI);
       //if (del.inc < 0.0) del.inc += 2.0*PI;
       //deltaOm = 2.0*PI*sqrt(delta/tRRv)/3.0;
       //del.lan += (2.0*drand48() - 1.0)*deltaOm;  /* lan random walks on tRRv */
       //del.lan = fmod(del.lan,2*PI);
       //if (del.lan < 0.0) del.lan += 2.0*PI;
       //deltaw = 2.0*PI*(delta/tprec)/3.0;
       //del.lop -= deltaw;  /* lop precesses due to mean stellar field */
       //if (del.lop < 0.0) del.lop += 2.0*PI;
       //del.lop = fmod(del.lop,2*PI);
       deltaM = 2.0*PI*delta/torb;  
       del.mea += deltaM;  /* update mean anomaly */
       del.mea = fmod(del.mea,2*PI);  /* Make sure between 0 and 2pi */
       fprintf(elem,"%g %g %g %g %g %g %g\n",t,del.sma,del.ecc,del.inc,del.lan,del.lop,del.mea);
       deltohel(G*(Mbh+Mbin),&del,&hel);  /* convert to cartesian */
       fprintf(orbit,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
	       t,hel.x,hel.y,hel.z,hel.vx,hel.vy,hel.vz);
     }

   fclose(orbit);
   fclose(elem);

}
