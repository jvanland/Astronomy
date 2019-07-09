#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10000

main()
{
   int i;
   double j0,j1,j2,j3;
   double x,y,z,vx,vy,vz,v;
   double j10,j11,j12,j13;
   FILE *initial;

   initial=fopen("initial.dat.bt","r");
   for (i=0; i<N; i++)
   {
     fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
       &j0,&j1,&j2,&j3,&x,&y,&z,&vx,&vy,&vz,&j10,&j11,&j12,&j13);
     v=sqrt(vx*vx+vy*vy+vz*vz);
     printf("%lg\n",v);
   }
   fclose(initial);
}
