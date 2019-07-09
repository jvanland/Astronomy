/* 
Draw a binary.
For use in hndrag.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "delaunay.h"
#include "helio.h"

#define G 1.0 /* Newton's constant for AU,yr,Msun*/
#define M 1.0e5    /* Mass of the more massive member of the binary */
#define m 1.0      /* Mass of the less massive member of the binary */
#define PI 3.141592653589793  /* pi to unnecessary precision */

void deltohel(double,struct delaunay *,struct helio *);

main()
{
  struct helio hel1,hel2;
  struct delaunay del;
  FILE *initial;
   
  initial=fopen("binary.dat","w");
  /* Define orbit of binary */
  del.sma=2000.0;
  del.ecc=0.99;
  del.inc=0.0;
  del.lan=0.4;
  del.lop=0.4;
  del.mea=3.0;
  /* Convert to cartesian coordinates and put c.o.m. at 0 */
  deltohel(G*(M+m),&del,&hel1);
  hel1.x*=m/(M+m);
  hel2.x=-(M/m)*hel1.x;
  hel1.y*=m/(M+m);
  hel2.y=-(M/m)*hel1.y;
  hel1.z*=m/(M+m);
  hel2.z=-(M/m)*hel1.z;
  hel1.vx*=m/(M+m);
  hel2.vx=-(M/m)*hel1.vx;
  hel1.vy*=m/(M+m);
  hel2.vy=-(M/m)*hel1.vy;
  hel1.vz*=m/(M+m);
  hel2.vz=-(M/m)*hel1.vz;

  fprintf(initial,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
	  M,hel1.x,hel1.y,hel1.z,hel1.vx,hel1.vy,hel1.vz);
  fprintf(initial,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
	  m,hel2.x,hel2.y,hel2.z,hel2.vx,hel2.vy,hel2.vz);

  fclose(initial);
}

