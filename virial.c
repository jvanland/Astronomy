/* Code to check whether the positions and velocities we
have generated come close to satisfying the virial theorem,
which says that for systems in dynamic equilibrium the
total kinetic energy equals -0.5 times the total potential
energy. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 1.0
#define PI 3.141592653589793
#define N 1001

main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Initfile.dat\n", argv[0]);
    return 1;
  }
   char *Initfile = argv[1];
   int i,j;
   double m[N],x[N]={0.0},y[N]={0.0},z[N]={0.0},d;
   double vx[N]={0.0},vy[N]={0.0},vz[N]={0.0};
   double kinetic, potential;
   FILE *initial;

   initial=fopen(Initfile,"r");

   kinetic=0.0;
   for (i=0; i<N; i++)
   {
      fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg",
         &m[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i]);
      kinetic+=0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
   }
   fclose(initial);
   potential=0.0;
   for (i=0; i<N-1; i++)
   {
      for (j=i+1; j<N; j++)
      {
         d=sqrt((x[j]-x[i])*(x[j]-x[i])+(y[j]-y[i])*(y[j]-y[i])+(z[j]-z[i])*(z[j]-z[i]));
         potential-=m[i]*m[j]/d;
      }
   }
   printf("kinetic=%lg, potential=%lg, T+U/2=%lg\n",
      kinetic,potential,kinetic+0.5*potential);
}
