#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10

main(){

   int i,j;
   double j0,j1,j2;
   double m[N],x[N]={0.0},y[N]={0.0},z[N]={0.0},d;
   double vx[N]={0.0},vy[N]={0.0},vz[N]={0.0};
   double j10,j11,j12,j13;
   double kinetic;
   FILE *initial;

   initial=fopen("ss.10000.bt","r");

   kinetic=0.0;
   for (i=0; i<N; i++)
   {
      fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
         &j0,&j1,&m[i],&j2,&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i],&j10,&j11,&j12,&j13);

      kinetic+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
      printf("%lg\n", kinetic);
   }
   fclose(initial);
}
