#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793

main()
{

  int i,t;
  double m, x, y, z, vx, vy, vz; 
  FILE *data;

  data=fopen("check.dat","r");

  for (i=0; i<10; i++)
    {
      fscanf(data,"%d %lg %lg %lg %lg %lg %lg %lg",
	     &t,&m,&x,&y,&z,&vx,&vy,&vz);

      vx *= 2.0*PI;
      vy *= 2.0*PI;
      vz *= 2.0*PI;

      printf("%g %.16g %.16g %.16g %.16g %.16g %.16g\n",m,x,y,z,vx,vy,vz);
    }
  fclose(data);
}
