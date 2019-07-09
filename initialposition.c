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
   FILE *dbin,*dsuper,*numbers;

   dbin=fopen("initial0.dat","r");
   dsuper=fopen(superdata,"r");

   fscanf(dbin,"%lg %lg %lg %lg %lg %lg %lg",
      &m1,&x1,&y1,&z1,&vx1,&vy1,&vz1);
   fscanf(dbin,"%lg %lg %lg %lg %lg %lg %lg",
      &m2,&x2,&y2,&z2,&vx2,&vy2,&vz2);
   xcm=(m1*x1+m2*x2)/(m1+m2);
   ycm=(m1*y1+m2*y2)/(m1+m2);
   zcm=(m1*z1+m2*z2)/(m1+m2);
   vxcm=(m1*vx1+m2*vx2)/(m1+m2);
   vycm=(m1*vy1+m2*vy2)/(m1+m2);
   vzcm=(m1*vz1+m2*vz2)/(m1+m2);
   printf("%lg %lg %lg %lg %lg %lg\n",xcm,ycm,zcm,vxcm,vycm,vzcm);
   fscanf(dsuper,"%d %lg %lg %lg %lg %lg %lg",
      &t,&xsuper,&ysuper,&zsuper,&vxsuper,&vysuper,&vzsuper);
   x1+=xsuper-xcm;
   x2+=xsuper-xcm;
   vx1+=vxsuper-vxcm;
   vx2+=vxsuper-vxcm;
   y1+=ysuper-ycm;
   y2+=ysuper-ycm;
   vy1+=vysuper-vycm;
   vy2+=vysuper-vycm;
   z1+=zsuper-zcm;
   z2+=zsuper-zcm;
   vz1+=vzsuper-vzcm;
   vz2+=vzsuper-vzcm;
   vx1*=2.0*PI;
   vx2*=2.0*PI;
   vy1*=2.0*PI;
   vy2*=2.0*PI;
   vz1*=2.0*PI;
   vz2*=2.0*PI;

   numbers=fopen("numbers.hnb","w");
   fprintf(numbers,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
      m1,x1,y1,z1,vx1,vy1,vz1);
   fprintf(numbers,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
      m2,x2,y2,z2,vx2,vy2,vz2);
   fclose(numbers);
   fclose(dbin);
   fclose(dsuper);
}
