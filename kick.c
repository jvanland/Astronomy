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
  double Minternal,dist,vmag1,vmag2;
  double x,y,z,vx,vy,vz;
  double j0,j1,j2,j3,j10,j11,j12,j13;
  FILE *init;
  char cmd[50];

  N_steps = 5000;

  //  x = 85.463919;
  //  y = -10.332869;
  //  z = 84.737564;
  //  dist = x*x + y*y + z*z;
  //  dist = sqrt(dist);
  //  vx = -2.697837;
  //  vy = 0.291634;
  //  vz = -1.935681;
  x = 9000;
  y = 0;
  z = 0;
  vx = 0;
  vy = 0;
  vz = 0;

  Minternal = 1.0e3 + 1.0e3*dist/(9e3);

  system("rm ss.*");
  init=fopen("initial.bt","w");
  fprintf(init,"%i %i %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %i\n",
  	  0,0,1.0,0.01,x,y,z,vx,vy,vz,0.0,0.0,0.0,3);
  fclose(init);
  system("cp initial.bt output2.dat");
  system("bt2ss initial.bt");

  sprintf(cmd,"cat template | sed 's/DCENTMASS/'%.16g'/' > ss.par", Minternal);
  system(cmd);
  system("./pkdgrav ss.par >& output2.log");

  for (i=1; i<N_steps; i++)
    {
      if(i % 100 == 0) printf("%d\n",i);
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
      if (dist < 3.0) return 0;
      Minternal = 1.0e3 + 1.0e3*dist/(9.0e3);
      if (Minternal > 2.0e3) Minternal = 2.0e3;
      //  if (i == 1209) {
      //	vmag1 = vx*vx + vy*vy + vz*vz;
      //	vmag1 = sqrt(vmag1); 
      //	vmag2 = vmag1 + 0.25; //This is the kick
      //	vx *= vmag2/vmag1;
      //	vy *= vmag2/vmag1;
      //	vz *= vmag2/vmag1;
      //     	init=fopen("initial.bt","w");
      //	fprintf(init,"%i %i %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %i\n",
      //	0,0,1.0,0.01,x,y,z,vx,vy,vz,0.0,0.0,0.0,3);
      //	fclose(init);
      //	system("bt2ss initial.bt");
      //     }
      sprintf(cmd,"cat template | sed 's/DCENTMASS/'%.16g'/' > ss.par", Minternal);
      system(cmd);
      system("./pkdgrav ss.par >> output2.log 2>&1");
    }
}
