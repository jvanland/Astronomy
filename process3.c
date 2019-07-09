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
#define head 0  /* Number of lines in header */

void heltodel(double,struct helio *,struct delaunay *);

main()
{
   int i,N=0;
   double m1,m2,mtot,dt,t=0;
   struct helio hel1,hel2;
   struct delaunay del1,del2;
   double x1,y1,z1,vx1,vy1,vz1;
   double x2,y2,z2,vx2,vy2,vz2;
   double mutinc,Lxbin,Lybin,Lzbin,Lbin;
   double Lxcm,Lycm,Lzcm,Lcm;
   FILE *data1,*data2;
   FILE *out1,*out2;
   char buf[5000],line[80];

   out1=fopen("comelem.dat","w");
   out2=fopen("binelem.dat","w");
   data1=fopen("fullbody1.dat","r");
   data2=fopen("fullbody2.dat","r");
   while((fgets(buf,5000,data1)) != NULL) {
     N++;
   }
   rewind(data1);
   N-=head;
   /* Read headers */
   for (i=0; i<head;i++) {
     fgets (line,sizeof line,data1); 
     fgets (line,sizeof line,data2);
   }

   m1=m2=10;

   for (i=0; i<N; i++)
   {
      fscanf(data1,"%lg %lg %lg %lg %lg %lg %lg",
         &dt,&x1,&y1,&z1,&vx1,&vy1,&vz1);
      fscanf(data2,"%lg %lg %lg %lg %lg %lg %lg",
         &dt,&x2,&y2,&z2,&vx2,&vy2,&vz2);
      t+=dt;
      /* C.O.M. Orbit */
      hel1.x=(m1*x1+m2*x2)/(m1+m2);
      hel1.y=(m1*y1+m2*y2)/(m1+m2);
      hel1.z=(m1*z1+m2*z2)/(m1+m2);
      hel1.vx=(m1*vx1+m2*vx2)/(m1+m2);
      hel1.vy=(m1*vy1+m2*vy2)/(m1+m2);
      hel1.vz=(m1*vz1+m2*vz2)/(m1+m2);
      mtot=m1+m2+1.0e5;
      heltodel(G*mtot,&hel1,&del1);
      del1.inc = fmod(del1.inc,2*PI);
      del1.lan = fmod(del1.lan,2*PI);
      del1.lop = fmod(del1.lop,2*PI);
      fprintf(out1,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",t,del1.sma,del1.ecc,del1.inc,del1.lan,del1.lop,del1.mea);
      /* Binary Orbit */
      hel2.x=x1-x2;
      hel2.y=y1-y2;
      hel2.z=z1-z2;
      hel2.vx=vx1-vx2;
      hel2.vy=vy1-vy2;
      hel2.vz=vz1-vz2;
      mtot=m1+m2;
      heltodel(G*mtot,&hel2,&del2);
      del2.lan = fmod(del2.lan,2*PI);
      del2.lop = fmod(del2.lop,2*PI);
      /* Calculate mutual inclination */
      Lxbin=hel2.y*hel2.vz-hel2.z*hel2.vy;
      Lybin=hel2.z*hel2.vx-hel2.x*hel2.vz;
      Lzbin=hel2.x*hel2.vy-hel2.y*hel2.vx;
      Lbin=sqrt(Lxbin*Lxbin+Lybin*Lybin+Lzbin*Lzbin);
      Lxbin/=Lbin;
      Lybin/=Lbin;
      Lzbin/=Lbin;
      Lxcm=hel1.y*hel1.vz-hel1.z*hel1.vy;
      Lycm=hel1.z*hel1.vx-hel1.x*hel1.vz;
      Lzcm=hel1.x*hel1.vy-hel1.y*hel1.vx;
      Lcm=sqrt(Lxcm*Lxcm+Lycm*Lycm+Lzcm*Lzcm);
      Lxcm/=Lcm;
      Lycm/=Lcm;
      Lzcm/=Lcm;
      mutinc=acos(Lxbin*Lxcm+Lybin*Lycm+Lzbin*Lzcm);
      fprintf(out2,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",t,del2.sma,del2.ecc,mutinc,del2.lan,del2.lop,del2.mea);
   }
   fclose(data1);
   fclose(data2);
}
