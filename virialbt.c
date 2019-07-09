/* Code to check whether the positions and velocities we
have generated come close to satisfying the virial theorem,
which says that for systems in dynamic equilibrium the
total kinetic energy equals -0.5 times the total potential
energy. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.673e-8
#define Msun 1.989e33
#define MSMBH 1.0e4   /* Mass of SMBH in solar masses */
#define AU 1.496e13
#define PI 3.141592653589793
#define yr2p (0.5*3.1557e7/PI)
#define N 1000

main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Initfile.bt\n", argv[0]);
    return 1;
  }
   char *Initfile = argv[1];
   int i,j;
   double j0,j1,j2;
   double m[N],x[N]={0.0},y[N]={0.0},z[N]={0.0},d;
   double vx[N]={0.0},vy[N]={0.0},vz[N]={0.0};
   double j10,j11,j12,j13;
   double kinetic, potential;
   FILE *initial;

   initial=fopen(Initfile,"r");

   kinetic=0.0;
   for (i=0; i<N; i++)
   {
      fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
         &j0,&j1,&m[i],&j2,&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i],&j10,&j11,&j12,&j13);
      m[i]*=Msun;
      x[i]*=AU;
      y[i]*=AU;
      z[i]*=AU;
      vx[i]*=AU/yr2p;
      vy[i]*=AU/yr2p;
      vz[i]*=AU/yr2p;
      kinetic+=0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
   }
   fclose(initial);
   potential=0.0;
   for (i=0; i<N-1; i++)
   {
      for (j=i+1; j<N; j++)
      {
         d=sqrt((x[j]-x[i])*(x[j]-x[i])+(y[j]-y[i])*(y[j]-y[i])+(z[j]-z[i])*(z[j]-z[i]));
         potential-=G*m[i]*m[j]/d;
      }
   }
   for (j=0; j<N; j++)
   {
      d=sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]);
      potential-=G*m[j]*MSMBH*Msun/d;
   }
   printf("kinetic=%lg, potential=%lg, T+U/2=%lg\n",
      kinetic,potential,kinetic+0.5*potential);
}
