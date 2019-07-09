/* Take a list of positions and velocities for the two
10 Msun black holes.  Using this, compute a and e
for the motion of the center of mass around the 10^5 Msun
black holes, a and e for the binary, and the relative
inclination of the binary to its superorbital plane. Units 
are AU, yr, and Msun*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "delaunay.h"
#include "helio.h"

#define G 39.47841760435743  /* Gravitational constant for Msun, AU, yr */
#define PI 3.141592653589793  /* pi to unnecessary precision */

void heltodel(double,struct helio *,struct delaunay *);

main()
{
   int i,N=0;
   double mbin,mtot,step,t=0;
   struct helio hel;
   struct delaunay del;
   FILE *data,*out;
   char buf[5000];

   out=fopen("fullelem.dat","w");
   data=fopen("f9980.dat","r");
   while((fgets(buf,5000,data)) != NULL) {
     N++;
   }
   rewind(data);

   mbin=20;

   for (i=0; i<N; i++)
   {
      fscanf(data,"%lg %lg %lg %lg %lg %lg %lg",
         &step,&hel.x,&hel.y,&hel.z,&hel.vx,&hel.vy,&hel.vz);
      hel.vx*=2.0*PI;
      hel.vy*=2.0*PI;
      hel.vz*=2.0*PI;
      t=step*0.1987/(2.0*PI);
      mtot=mbin+1.0e5;
      heltodel(G*mtot,&hel,&del);
      del.inc = fmod(del.inc,2*PI);
      del.lan = fmod(del.lan,2*PI);
      del.lop = fmod(del.lop,2*PI);
      del.mea = fmod(del.mea,2*PI);
      fprintf(out,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",t,del.sma,del.ecc,del.inc,del.lan,del.lop,del.mea);
   }
   fclose(data);
   fclose(out);
}
