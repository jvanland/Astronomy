#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 6.673e-8      /* Newton's constant in cgs */
#define Msun 1.989e33   /* Solar mass in grams */
#define AU 1.496e13  /* AU in cm */
#define PI 3.141592653589793  /* pi to unnecessary precision */
#define yr2p (0.5*3.1557e7/PI)  /* One year over 2 pi, in seconds */
#define N_files 1099
#define N_obj 10000
#define N_bin 20

typedef struct
{
   double x;     /* x component of vector. */
   double y;     /* y component of vector. */
   double z;     /* z component of vector. */
} vector;

void zerovec(vector *r);    /* Set components of vector to zero */
vector vecadd(vector a, vector b);   /* Vector addition a+b. */
vector vecsub(vector a, vector b); /* Vector subtrction a-b */
vector vecmult(double d, vector a);  /* Multiplication: da. */
vector vecdiv(vector a, double d);  /* Division: a/d. */
vector crossprod(vector a, vector b); /* Cross product, axb */
double magnitude(vector a);  /* Magnitude of vector a */

main()
{

  int i,j,n;
  vector *r, *v;
  vector rdiff,vdiff,r_cm,a1,a2,a3;
  float h1,h2,h4,h11,h12,h13,h14;
  double M[N_obj],SuperM,Mtot,t;
  double magr,magv,magr_cm,Fbin,Ftide;
  double l,e,a4,a5,a6,r_p;
  char fp[50];
  FILE *output,*pkd;


  SuperM = 1.0e5*Msun;

  r=calloc(N_obj,sizeof(vector));
  v=calloc(N_obj,sizeof(vector));

   
  output=fopen("close.dat","w");

  for (i=0; i<N_files+1; i++)
    {
      if (3*i<10) sprintf(fp,"ss.00000%d000.bt",3*i);
      else if (3*i<100) sprintf(fp,"ss.0000%d000.bt",3*i);
      else if (3*i<1000) sprintf(fp,"ss.000%d000.bt",3*i);
      else sprintf(fp,"ss.00%d000.bt",3*i);
      pkd=fopen(fp,"r");

      printf("%d\n",3000*i);

      t=i*3000*0.1987/(2.0*PI);

      for (j=0; j<N_obj; j++)
	{
	  fscanf(pkd,"%f %f %lg %f %lg %lg %lg %lg %lg %lg %f %f %f %f",
		 &h1,&h2,&M[j],&h4,&r[j].x,&r[j].y,&r[j].z,&v[j].x,&v[j].y,&v[j].z,&h11,&h12,&h13,&h14);
	}
      fclose(pkd);

      for (n=0; n<N_bin; n++)
	{
	  for (j=0; j<N_obj; j++)
	    {
	      if (j != 9980+n) {
		Mtot = (M[9980+n]+M[j])*Msun;
		rdiff = vecsub(r[9980+n],r[j]);
		magr = magnitude(rdiff)*AU;
		a1 = vecmult(M[9980+n]*Msun,r[9980+n]);
		a2 = vecmult(M[j]*Msun,r[j]);
		a3 = vecadd(a1,a2);
		r_cm = vecdiv(a3,Mtot);
		magr_cm = magnitude(r_cm)*AU;
		Fbin = G*Mtot/(magr*magr);
		Ftide = 2*G*SuperM*magr/(magr_cm*magr_cm*magr_cm);

		if (Ftide < Fbin) {
		  vdiff = vecsub(v[9980+n],v[j]);
		  magv = magnitude(vdiff)*AU/yr2p;
		  l = magnitude(crossprod(rdiff,vdiff))*AU*AU/yr2p;
		  e = 0.5*magv*magv - G*Mtot/magr;
		  a4 = G*Mtot/e;
		  a5 = 2*l*l/e;
		  a6 = sqrt(a4*a4+a5);
		  r_p = (-a4 + a6)/(2*AU);
		  fprintf(output,"%d\t%d\t%lg\t%lg\t%lg\t%lg\n",9980+n,j,Fbin/Ftide,r_p,magr/AU,t);
	      }
	      }
		}
	}
    }

  fclose(output);
  free(r);
  free(v);
}

void zerovec(vector *r)
{
   (*r).x=0.0;
   (*r).y=0.0;
   (*r).z=0.0;
}

vector vecadd(vector a, vector b)
{
   vector sum;

   sum.x=a.x+b.x;
   sum.y=a.y+b.y;
   sum.z=a.z+b.z;
   return sum;
}

vector vecsub(vector a, vector b)
{
   vector diff;

   diff.x=a.x-b.x;
   diff.y=a.y-b.y;
   diff.z=a.z-b.z;
   return diff;
}

vector vecmult(double d, vector a)
{
   vector mult;

   mult.x=d*a.x;
   mult.y=d*a.y;
   mult.z=d*a.z;
   return mult;
}

vector vecdiv(vector a, double d)
{
   vector div;

   div.x=a.x/d;
   div.y=a.y/d;
   div.z=a.z/d;
   return div;
}

vector crossprod(vector a, vector b)
{
   vector c;

   c.x=a.y*b.z-a.z*b.y;
   c.y=a.z*b.x-a.x*b.z;
   c.z=a.x*b.y-a.y*b.x;
   return c;
}

double magnitude(vector a)
{
   double mag=0.0;

   mag=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
   return mag;
}
