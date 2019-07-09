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
   double kinetic;
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
      kinetic=0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
      printf("%d %lg\n",i,kinetic);
   }
   fclose(initial);
   
}
