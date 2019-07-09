/* Code to run a series of 1 step , 1 particle 
pkdgrav runs where the mass of the SMBH changes 
depending on the distance of the particle. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 1.0  /* Newton's constant in pkdgrav units */
#define PI 3.141592653589793  /* pi to unnecessary precision */

main()
{
  int i,N_steps;
  double a,e,Minternal,dist;
  double x,y,z,vx,vy,vz;
  double j0,j1,j2,j3,j10,j11,j12,j13;
  FILE *init;
  char cmd[50];

  a = 6000.0;
  //  e = 0.998;
  N_steps = 3000;

  x = 100.0131091;
  Minternal = 1.0e5 + 1.0e4*x/(1.998*a);
  //  y = 0.0;
  //  z = 0.0;
  //  vx = 0.0;
  //  vy = G*Minternal*((2.0/x) - (1.0/a));
  //  vy = sqrt(vy);
  //  vz = 0.0;

  //  system("rm ss.*");
  //  init=fopen("initial.bt","w");
  //  fprintf(init,"%i %i %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %i\n",
  //	  0,0,1.0,1.0,x,y,z,vx,vy,vz,0.0,0.0,0.0,3);
  //  fclose(init);
  system("cp initial.bt output2.dat");
  system("bt2ss initial.bt");

  sprintf(cmd,"cat template | sed 's/DCENTMASS/'%.16g'/' > ss.par", Minternal);
  system(cmd);
  system("./pkdgrav ss.par >& output2.log");

  for (i=1; i<N_steps; i++)
    {
      system("mv ss.1 initial.ss");
      system("rm ss.*");
      system("ss2bt initial.ss");
      system("cat initial.bt >> output2.dat");
      init=fopen("initial.bt","r");
      fscanf(init,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
	     &j0,&j1,&j2,&j3,&x,&y,&z,&vx,&vy,&vz,&j10,&j11,&j12,&j13);
      fclose(init);
      dist = x*x + y*y + z*z;
      dist = sqrt(dist);
      Minternal = 1.0e5 + 1.0e4*dist/(1.998*a);
      sprintf(cmd,"cat template | sed 's/DCENTMASS/'%.16g'/' > ss.par", Minternal);
      system(cmd);
      system("./pkdgrav ss.par >> output2.log");
    }
}
