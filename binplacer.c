/* Take the position and velocity of two objects relative
to their center of mass, then place this on top of the
motion of the center of mass. 
Called by driver.c
*/

#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265359

main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Data_file\n", argv[0]);
    return 1;
  }
   char *superdata = argv[1];

   int t;
   double m1,x1,y1,z1,vx1,vy1,vz1;
   double m2,x2,y2,z2,vx2,vy2,vz2;
   double xcm,ycm,zcm,vxcm,vycm,vzcm;
   double xsuper,ysuper,zsuper,vxsuper,vysuper,vzsuper;
   FILE *dbin,*dsuper,*hnb;

   dbin=fopen("binary.dat","r");
   dsuper=fopen(superdata,"r");

   fscanf(dbin,"%lg %lg %lg %lg %lg %lg %lg",
      &m1,&x1,&y1,&z1,&vx1,&vy1,&vz1);
   fscanf(dbin,"%lg %lg %lg %lg %lg %lg %lg",
      &m2,&x2,&y2,&z2,&vx2,&vy2,&vz2);
   fscanf(dsuper,"%d %lg %lg %lg %lg %lg %lg",
      &t,&xsuper,&ysuper,&zsuper,&vxsuper,&vysuper,&vzsuper);
   x1+=xsuper;
   x2+=xsuper;
   vx1+=vxsuper;
   vx2+=vxsuper;
   y1+=ysuper;
   y2+=ysuper;
   vy1+=vysuper;
   vy2+=vysuper;
   z1+=zsuper;
   z2+=zsuper;
   vz1+=vzsuper;
   vz2+=vzsuper;

   hnb=fopen("binary.hnb","w");
   fprintf(hnb,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
      m1,x1,y1,z1,vx1,vy1,vz1);
   fprintf(hnb,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
      m2,x2,y2,z2,vx2,vy2,vz2);
   fclose(hnb);
   fclose(dbin);
   fclose(dsuper);
}
