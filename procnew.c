/* Take a list of positions and velocities for the two
10 Msun black holes.  Using this, compute a and e
for the motion of the center of mass around the 10^4 Msun
black holes, a and e for the binary, and the relative
inclination of the binary to its superorbital plane. Units 
are AU, yr, and Msun*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "helio.h"
#include "delaunay.h"

void heltodel(double mu,struct helio *h,struct delaunay *d);

main()
{
  int i,N=0;
   double m1,m2,mu,t;
   double x1,y1,z1,vx1,vy1,vz1;
   double x2,y2,z2,vx2,vy2,vz2;
   struct helio bCoords, cCoords;
   struct delaunay bElem, cElem;
   FILE *data1,*data2;
   char buf[5000];

   data1=fopen("body1.dat","r");
   data2=fopen("body2.dat","r");
   while((fgets(buf,5000,data1)) != NULL) {
     N++;
   }
   rewind(data1);

   m1=m2=10;

   for (i=0; i<N; i++)
   {
      fscanf(data1,"%lg %lg %lg %lg %lg %lg %lg",
         &t,&x1,&y1,&z1,&vx1,&vy1,&vz1);
      fscanf(data2,"%lg %lg %lg %lg %lg %lg %lg",
         &t,&x2,&y2,&z2,&vx2,&vy2,&vz2);
      cCoords.x=(m1*x1+m2*x2)/(m1+m2);
      cCoords.y=(m1*y1+m2*y2)/(m1+m2);
      cCoords.z=(m1*z1+m2*z2)/(m1+m2);
      cCoords.vx=(m1*vx1+m2*vx2)/(m1+m2);
      cCoords.vy=(m1*vy1+m2*vy2)/(m1+m2);
      cCoords.vz=(m1*vz1+m2*vz2)/(m1+m2);
      bCoords.x=x1-x2;
      bCoords.y=y1-y2;
      bCoords.z=z1-z2;
      bCoords.vx=vx1-vx2;
      bCoords.vy=vy1-vy2;
      bCoords.vz=vz1-vz2;
      mu=(4*CONST_PI*CONST_PI)*(m1+m2);  // G*Mtot in these units
      heltodel(mu,&bCoords,&bElem);
      mu=(4*CONST_PI*CONST_PI)*(m1+m2+1.0e5);
      heltodel(mu,&cCoords,&cElem);
      printf("%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
	     t,bElem.sma,bElem.ecc,bElem.sma*(1.0-bElem.ecc),cElem.sma,cElem.ecc,bElem.inc-cElem.inc);
   }
   fclose(data1);
   fclose(data2);
}

