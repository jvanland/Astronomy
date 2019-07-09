/* Sets up an orbit around a SMBH where
one or more of the orbital elements changes
predicatbly over time.  Prints out the
cartesian coordinates of the orbit.  For
use with hnbody to examine Kozai. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "delaunay.h"
#include "helio.h"

#define G 39.47841760435743  /* Gravitational constant for Msun, AU, yr */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define Mstar 20.0  /* Mass of star */
#define Mbh 1.0e5  /* Mass of SMBH */
#define delta 10.0  /* Time between outputs (yr) */
#define Tfin 1.0e5  /* Total time */


void deltohel(double,struct delaunay *,struct helio *);

main(int argc, char *argv[ ])
{
  if(argc != 4) {
    printf("Usage: %s kep_elem(1-5) init_val grow/decay(1,2)\n", argv[0]);
    return 1;
  }

   int i,N;
   struct helio hel;
   struct delaunay del;
   double deltaM, t=0;
   int elem = atoi(argv[1]);
   double change = atof(argv[2]);
   int type = atoi(argv[3]);
   FILE *orbit,*init;

   orbit=fopen("orbit.dat","w");

   /*Read in initial orbital elements */
   init=fopen("eleminit.dat","r");
   fscanf(init,"%lg %lg %lg %lg %lg %lg",
      &del.sma,&del.ecc,&del.inc,&del.lan,&del.lop,&del.mea);
   fclose(init);
   if (elem == 1) del.sma = change;
   else if (elem == 2) del.ecc = change;
   else if (elem == 3) del.inc = change;
   else if (elem == 4) del.lan = change;
   else if (elem == 5) del.lop = change;
   /* Convert to cartesian coordinates */
   deltohel(G*(Mbh+Mstar),&del,&hel);
   fprintf(orbit,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
	   t,hel.x,hel.y,hel.z,hel.vx,hel.vy,hel.vz);

   N = Tfin/delta;

   for (i=0; i<N; i++)
     {
       t += delta;  /* update time */
       if (type == 1)  /* Linear Growth */
	 {
	   if (elem == 1) del.sma = change*2*(t/Tfin);
	   else if (elem == 2) del.ecc = change*2*(t/Tfin);
	   else if (elem == 3) del.inc = change+2*PI*(t/Tfin);
	   else if (elem == 4) del.lan = change+2*PI*(t/Tfin);
	   else if (elem == 5) del.lop = change+2*PI*(t/Tfin);
	 }
       if (type == 2)  /* Linear Decay */
	 {
	   if (elem == 1) del.sma = change*2*(-t/Tfin);
	   else if (elem == 2) del.ecc = change*2*(-t/Tfin);
	   else if (elem == 3) del.inc = change-2*PI*(t/Tfin);
	   else if (elem == 4) del.lan = change-2*PI*(t/Tfin);
	   else if (elem == 5) del.lop = change-2*PI*(t/Tfin);
	 }
       deltaM = sqrt((G*(Mbh+Mstar))/(del.sma*del.sma*del.sma))*delta;  
       del.mea += deltaM;  /* update mean anomaly */
       del.mea = fmod(del.mea,2*PI);  /* Make sure between 0 and 2pi */
       deltohel(G*(Mbh+Mstar),&del,&hel);  /* convert to cartesian */
       fprintf(orbit,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
	       t,hel.x,hel.y,hel.z,hel.vx,hel.vy,hel.vz);
     }

   fclose(orbit);

}
