/* Prints out the net momentum and angular momentum for a
collection of objects. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 1.0      /* Newton's constant in units of Msun, Au, yr/2pi */
#define Msun 1.0   /* Msun in Msun */
#define M 1.0e3   /* Mass of SMBH in solar masses */
#define N_obj 1000     /* Number of stars */
#define N_steps 100000    /* Number of pkd steps */
#define Out_int 100    /* Interval of pkd output in steps */
#define PI 3.141592653589793  /* pi to unnecessary precision */

typedef struct
{
   double x;     /* x component of vector. */
   double y;     /* y component of vector. */
   double z;     /* z component of vector. */
} vector;

void zerovec(vector *r);    /* Set components of vector to zero */
vector vecadd(vector a, vector b);   /* Vector addition a+b. */
vector vecmult(double d, vector a);  /* Multiplication: da. */
vector crossprod(vector a, vector b); /* Cross product, axb */
double magnitude(vector a);  /* Magnitude of vector a */

main(int argc, char *argv[ ])
{
//  if(argc != 2) {
//    printf("Usage: %s Initfile.bt\n", argv[0]);
//    return 1;
//  }

  int i,j,k;
  double t0,t1,t2;
  double m[N_obj],L;
  vector *r,*v;
  double t10,t11,t12,t13;
  vector mom,angmom;
  double kinetic, potential, tot, d;
  char fp[50];
  FILE *pkd;

  r=calloc(N_obj,sizeof(vector));
  v=calloc(N_obj,sizeof(vector));

  for (j=Out_int; j<=N_steps; j+=Out_int)
    {
      if (j<1000) sprintf(fp,"ss.000%d.bt",j);
      else if (j<10000) sprintf(fp,"ss.00%d.bt",j);
      else if (j<100000) sprintf(fp,"ss.0%d.bt",j);
      else sprintf(fp,"ss.%d.bt",j);
      pkd=fopen(fp,"r");

      zerovec(&mom);
      zerovec(&angmom);
      kinetic=0.0;

      for (i=0; i<N_obj; i++)
	{
	  fscanf(pkd,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
		 &t0,&t1,&m[i],&t2,&r[i].x,&r[i].y,&r[i].z,&v[i].x,&v[i].y,&v[i].z,&t10,&t11,&t12,&t13);
	  mom=vecadd(mom,vecmult(m[i],v[i]));
	  angmom=vecadd(angmom,vecmult(m[i],crossprod(r[i],v[i])));
	  kinetic+=0.5*m[i]*(v[i].x*v[i].x+v[i].y*v[i].y+v[i].z*v[i].z);
	}
      L=magnitude(angmom);
      fclose(pkd);

      potential=0.0;
      for (i=0; i<N_obj-1; i++)
	{
	  for (k=i+1; k<N_obj; k++)
	    {
	      d=sqrt((r[k].x-r[i].x)*(r[k].x-r[i].x)+(r[k].y-r[i].y)*(r[k].y-r[i].y)+(r[k].z-r[i].z)*(r[k].z-r[i].z));
	      potential-=G*m[i]*m[k]/d;
	    }
	}
      for (k=0; k<N_obj; k++)
	{
	  d=sqrt(r[k].x*r[k].x+r[k].y*r[k].y+r[k].z*r[k].z);
	  potential-=G*m[k]*M*Msun/d;
	}
      tot = kinetic + potential;

      printf("%d %lg %lg %lg %lg %lg %lg\n",j,tot,kinetic+0.5*potential,L,mom.x,mom.y,mom.z);
    }

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

vector vecmult(double d, vector a)
{
   vector mult;

   mult.x=d*a.x;
   mult.y=d*a.y;
   mult.z=d*a.z;
   return mult;
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
