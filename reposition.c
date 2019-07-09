/* Take the position and velocity of two objects relative
to their center of mass, then place this on top of the
motion of the center of mass. 
Called by driver.c, driver2.c
*/

#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265359

main()
{
   int t;
   double m1,x1,y1,z1,vx1,vy1,vz1;
   double m2,x2,y2,z2,vx2,vy2,vz2;
   double xcm,ycm,zcm,vxcm,vycm,vzcm;
   double xsuper,ysuper,zsuper,vxsuper,vysuper,vzsuper;
   FILE *dbin,*dsuper,*numbers;

   dbin=fopen("initial.dat","r");
   dsuper=fopen("super.dat","r");

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
   fscanf(dsuper,"%d %lg %lg %lg %lg %lg %lg",
      &t,&xsuper,&ysuper,&zsuper,&vxsuper,&vysuper,&vzsuper);
   x1+=xsuper-xcm;
   x2+=xsuper-xcm;
   vx1+=vxsuper*2.0*PI-vxcm;
   vx2+=vxsuper*2.0*PI-vxcm;
   y1+=ysuper-ycm;
   y2+=ysuper-ycm;
   vy1+=vysuper*2.0*PI-vycm;
   vy2+=vysuper*2.0*PI-vycm;
   z1+=zsuper-zcm;
   z2+=zsuper-zcm;
   vz1+=vzsuper*2.0*PI-vzcm;
   vz2+=vzsuper*2.0*PI-vzcm;

   numbers=fopen("numbers.hnb","w");
   fprintf(numbers,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
      m1,x1,y1,z1,vx1,vy1,vz1);
   fprintf(numbers,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
      m2,x2,y2,z2,vx2,vy2,vz2);
   fclose(numbers);
   fclose(dbin);
   fclose(dsuper);
}
