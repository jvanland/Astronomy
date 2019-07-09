/* Take a list of positions and velocities for the two
10 Msun black holes.  Using this, compute a and e
for the motion of the center of mass around the 10^5 Msun
black holes, a and e for the binary, and the relative
inclination of the binary to its superorbital plane. Units 
are AU, yr, and Msun*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793
#define G 39.47841760435743
#define MSMBH 1.0e5

void aei(double Mtot, double x, double y, double z, double vx, 
double vy, double vz, double *a, double *e, double *Lx, double *Ly, 
double *Lz);

main()
{
   int i,N=0;
   double m1,m2,mtot,t;
   double x1,y1,z1,vx1,vy1,vz1;
   double x2,y2,z2,vx2,vy2,vz2;
   double xcm,ycm,zcm,vxcm,vycm,vzcm;
   double abin,ebin,Lxbin,Lybin,Lzbin;
   double acm,ecm,Lxcm,Lycm,Lzcm;
   double inc;
   FILE *data1,*data2;
   char buf[5000],line[80];

   data1=fopen("body1.dat","r");
   data2=fopen("body2.dat","r");
   while((fgets(buf,5000,data1)) != NULL) {
     N++;
   }
   rewind(data1);

   m1=m2=10;
   // Manually input N
   N = 3900;

   for (i=0; i<N; i++)
   {
      fscanf(data1,"%lg %lg %lg %lg %lg %lg %lg",
         &t,&x1,&y1,&z1,&vx1,&vy1,&vz1);
      fscanf(data2,"%lg %lg %lg %lg %lg %lg %lg",
         &t,&x2,&y2,&z2,&vx2,&vy2,&vz2);
      xcm=(m1*x1+m2*x2)/(m1+m2);
      ycm=(m1*y1+m2*y2)/(m1+m2);
      zcm=(m1*z1+m2*z2)/(m1+m2);
      vxcm=(m1*vx1+m2*vx2)/(m1+m2);
      vycm=(m1*vy1+m2*vy2)/(m1+m2);
      vzcm=(m1*vz1+m2*vz2)/(m1+m2);
      x1-=xcm;
      x2-=xcm;
      y1-=ycm;
      y2-=ycm;
      z1-=zcm;
      z2-=zcm;
      vx1-=vxcm;
      vx2-=vxcm;
      vy1-=vycm;
      vy2-=vycm;
      vz1-=vzcm;
      vz2-=vzcm;
      mtot=m1+m2;
      aei(mtot,x1-x2,y1-y2,z1-z2,vx1-vx2,vy1-vy2,vz1-vz2,&abin,&ebin,&Lxbin,&Lybin,&Lzbin);
      mtot=m1+m2+MSMBH;
      aei(mtot,xcm,ycm,zcm,vxcm,vycm,vzcm,&acm,&ecm,&Lxcm,&Lycm,&Lzcm);
      inc=acos(Lxbin*Lxcm+Lybin*Lycm+Lzbin*Lzcm);
      printf("%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",t,abin,ebin,abin*(1.0-ebin),acm,ecm,inc);
   }
   fclose(data1);
   fclose(data2);
}

void aei(double Mtot, double x, double y, double z, double vx, double vy, 
double vz, double *a, double *e, double *Lx, double *Ly, double *Lz)
{
/* Given position and velocity, compute semimajor axis, eccentricity,
and inclination */

   double dist,vel,energy,L;

   dist=sqrt(x*x+y*y+z*z);
   vel=sqrt(vx*vx+vy*vy+vz*vz);
   energy=0.5*vel*vel-G*Mtot/dist;
   *a=-G*Mtot/(2.0*energy);
   *Lx=y*vz-z*vy;
   *Ly=z*vx-x*vz;
   *Lz=x*vy-y*vx;
   L=sqrt((*Lx)*(*Lx)+(*Ly)*(*Ly)+(*Lz)*(*Lz));
   *Lx/=L;
   *Ly/=L;
   *Lz/=L;
   *e=sqrt(1.0-L*L/(G*Mtot*(*a)));
}
