#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define AU 1.496e13  /* AU in cm */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define yr2p (0.5*3.1557e7/PI)  /* One year over 2 pi, in seconds */
#define N_files 1100
#define N_obj 10000

typedef struct
{
   double x;     /* x component of vector. */
   double y;     /* y component of vector. */
   double z;     /* z component of vector. */
} vector;

double magnitude(vector a);  /* Magnitude of vector a */

main()
{

  int i,j;
  vector *r;
  float h1,h2,h3,h4,h8,h9,h10,h11,h12,h13,h14;
  double magr,t;
  char fp[50];
  FILE *output,*pkd;

  r=calloc(N_obj,sizeof(vector));
   
  output=fopen("radius.dat","w");

  for (i=1; i<+N_files; i++)
    {
      if (3*i<10) sprintf(fp,"ss.000%d0.bt",3*i);
      else if (3*i<100) sprintf(fp,"ss.00%d0.bt",3*i);
      else if (3*i<1000) sprintf(fp,"ss.0%d0.bt",3*i);
      else sprintf(fp,"ss.%d0.bt",3*i);
      pkd=fopen(fp,"r");

      printf("%s\n",fp);

      t=i*30*0.1987/(2.0*PI);

      for (j=0; j<N_obj; j++)
	{
	  fscanf(pkd,"%f %f %lg %f %lg %lg %lg %lg %lg %lg %f %f %f %f",
		 &h1,&h2,&h3,&h4,&r[j].x,&r[j].y,&r[j].z,&h8,&h9,&h10,&h11,&h12,&h13,&h14);
	}
      fclose(pkd);

      for (j=0; j<N_obj; j++)
	{
	    magr = magnitude(r[j]);

	  if (magr < 50) {
	    fprintf(output,"%d\t%lg\t%lg\n",j,magr,t);
	  }
	}
    }

  fclose(output);
  free(r);
}

double magnitude(vector a)
{
   double mag=0.0;

   mag=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
   return mag;
}
