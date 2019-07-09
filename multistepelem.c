/* Converts list of Cartesian Coordinates in bt
format to Orbital elements. Units are AU, yr, and Msun*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "delaunay.h"
#include "helio.h"

#define G 39.47841760435743  /* Gravitational constant for Msun, AU, yr */
#define PI 3.141592653589793  /* pi to unnecessary precision */
//  #define M 1.0e5  /* Mass of SMBH */

void heltodel(double,struct helio *,struct delaunay *);

main()
{
   int i,N=0;
   double m,M,mtot,energy,dist;
   double j0,j1,j3,j10,j11,j12,j13;
   struct helio hel;
   struct delaunay del;
   FILE *data;
   char buf[5000];

   data=fopen("output2.dat","r");
   while((fgets(buf,5000,data)) != NULL) {
     N++;
   }
   rewind(data);

   for (i=0; i<N; i++)
   {
      fscanf(data,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
	     &j0,&j1,&m,&j3,&hel.x,&hel.y,&hel.z,&hel.vx,&hel.vy,&hel.vz,&j10,&j11,&j12,&j13);
      hel.vx*=2.0*PI;
      hel.vy*=2.0*PI;
      hel.vz*=2.0*PI;
      dist = hel.x*hel.x + hel.y*hel.y + hel.z*hel.z;
      dist = sqrt(dist);
      M = 1.0e3 + 1.0e3*dist/(9e3);
      mtot=m+M;
      heltodel(G*mtot,&hel,&del);
      del.inc = fmod(del.inc,2*PI);
      del.lan = fmod(del.lan,2*PI);
      del.lop = fmod(del.lop,2*PI);
      del.mea = fmod(del.mea,2*PI);
      energy = -G*M/(2.0*del.sma);
      printf("%.16g %.16g %.16g %.16g %.16g\n",dist,del.sma,del.ecc,del.inc,energy);
   }
   fclose(data);
}
