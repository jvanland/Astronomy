#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "delaunay.h"
#include "helio.h"

#define G 39.47841760435743  /* Gravitational constant for Msun, AU, yr */
#define PI 3.141592653589793  /* pi to unnecessary precision */

void fg(double,double *,double *,double);

main()
{
  double mu,x[3],v[3], dt, converge, dx;

  mu = 100001;
  x[0] = 0.3063133432227427;
  x[1] = 0.1029475296029963;
  x[2] = 0;
  v[0] = -714.1461849953769;
  v[1] = -94.38030391886426;
  v[2] = 0;
  dt = 6.467798563298969e-07;

  converge = 2.4543692606169356e-20;
  dx = 6.6733460567803446e-20;
  if (fabs(dx) <= converge) printf("yes\n");
  else printf("no\n"); 

  fg(mu,x,v,dt);

  printf("%.16g %.16g %.16g %.16g %.16g %.16g\n",
	 x[0],x[1],x[2],v[0],v[1],v[2]);


}
