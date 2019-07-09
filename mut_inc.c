/* Take a list of positions and velocities from
a 3-body hndrag run and convert to orbital elements 
including mutual inclination.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793
#define G 39.47841760435743 /* Newton's constant for AU,yr,Msun*/
#define M1 1.0 /* Mass of central object */
#define M2 1.0e-3
#define M3 1.0e-3

void aei(double Mtot, double x, double y, double z, double vx, 
double vy, double vz, double *a, double *e, double *Lx, double *Ly, 
double *Lz);

int main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Directory\n", argv[0]);
    return 1;
  }

   char *dir = argv[1];
   int i,N=0;
   double t;
   double x1,y1,z1,vx1,vy1,vz1;
   double x2,y2,z2,vx2,vy2,vz2;
   double x3,y3,z3,vx3,vy3,vz3;
   double ain,ein,Lxin,Lyin,Lzin;
   double aout,eout,Lxout,Lyout,Lzout;
   double inc_in,inc_out,inc_mut;
   FILE *data1,*data2,*data3;
   char buf[5000];

   char file1[50];
   sprintf(file1,"%s/body0.dat", dir);
   data1=fopen(file1,"r");
   char file2[50];
   sprintf(file2,"%s/body1.dat", dir);
   data2=fopen(file2,"r");
   char file3[50];
   sprintf(file3,"%s/body2.dat", dir);
   data3=fopen(file3,"r");

   while((fgets(buf,5000,data1)) != NULL) {
     N++;
   }
   rewind(data1);

   for (i=0; i<N; i++)
   {
      fscanf(data1,"%lg %lg %lg %lg %lg %lg %lg",
         &t,&x1,&y1,&z1,&vx1,&vy1,&vz1);
      fscanf(data2,"%lg %lg %lg %lg %lg %lg %lg",
         &t,&x2,&y2,&z2,&vx2,&vy2,&vz2);
      fscanf(data3,"%lg %lg %lg %lg %lg %lg %lg",
         &t,&x3,&y3,&z3,&vx3,&vy3,&vz3);

      aei(M1+M2,x2-x1,y2-y1,z2-z1,vx2-vx1,vy2-vy1,vz2-vz1,&ain,&ein,&Lxin,&Lyin,&Lzin);
      inc_in=acos(Lzin)*180.0/PI;
      aei(M1+M3,x3-x1,y3-y1,z3-z1,vx3-vx1,vy3-vy1,vz3-vz1,&aout,&eout,&Lxout,&Lyout,&Lzout);
      inc_out=acos(Lzout)*180.0/PI;
      inc_mut=acos(Lxin*Lxout+Lyin*Lyout+Lzin*Lzout)*180.0/PI;

      printf("%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",t,ain,aout,ein,eout,inc_in,inc_out,inc_mut);
   }
   fclose(data1);
   fclose(data2);
   fclose(data3);
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
