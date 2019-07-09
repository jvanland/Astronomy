/* Calculate net force on object near pericenter
at two timesteps. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 39.47841760435743  /* Gravitational constant for Msun, AU, yr */
#define PI 3.141592653589793  /* pi to unnecessary precision */

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
double dotprod(vector a, vector b);  /* Dot product a.b */
double magnitude(vector a);  /* Magnitude of vector a */

main()
{

  int i,N=10000;
  float h1,h2,h4,h11,h12,h13,h14;
  double delta,M[N],Mbh=1.0e5;
  double magr,factor,diffen,orben,ratio;
  vector force1,force2,rdiff;
  vector *r, *v;
  vector startr,endr,startv,endv;
  FILE *init1,*init2;

  r=calloc(N,sizeof(vector));
  v=calloc(N,sizeof(vector));

  delta = 0.031624087;

  startr.x = -9.210517834850513;
  startr.y = -3.894144481641693;
  startr.z = -0;
  startv.x = 340.3150794113257;
  startv.y = -804.9208557010577;
  startv.z = -154.5141356773306;
  endr.x = 10.88856056323082;
  endr.y = -19.26191437670382;
  endr.z = -3.886536865275175;
  endv.x = 557.659169958878;
  endv.y = -183.9186351506658;
  endv.z = -68.34766777386952;

  init1=fopen("start.bt","r");
  zerovec(&force1);

  for (i=0; i<N; i++)
    {
      fscanf(init1,"%f %f %lg %f %lg %lg %lg %lg %lg %lg %f %f %f %f",
	     &h1,&h2,&M[i],&h4,&r[i].x,&r[i].y,&r[i].z,&v[i].x,&v[i].y,&v[i].z,&h11,&h12,&h13,&h14);

      rdiff = vecsub(startr,r[i]);
      magr = magnitude(rdiff);
      factor = -G*M[i]/(magr*magr*magr);
      force1 = vecadd(force1,vecmult(factor,rdiff));
    }

  init2=fopen("end.bt","r");
  zerovec(&force2);

  for (i=0; i<N; i++)
    {
      fscanf(init2,"%f %f %lg %f %lg %lg %lg %lg %lg %lg %f %f %f %f",
	     &h1,&h2,&M[i],&h4,&r[i].x,&r[i].y,&r[i].z,&v[i].x,&v[i].y,&v[i].z,&h11,&h12,&h13,&h14);

      rdiff = vecsub(endr,r[i]);
      magr = magnitude(rdiff);
      factor = -G*M[i]/(magr*magr*magr);
      force2 = vecadd(force2,vecmult(factor,rdiff));
    }

  diffen = (dotprod(force1,startv) - dotprod(force2,endv))*delta;
  orben = G*Mbh/(2.0*2000.0);
  ratio = diffen/orben;

  printf("%g %g %g\n",force1.x,force1.y,force1.z);
  printf("%g %g %g\n",force2.x,force2.y,force2.z);
  printf("%g %g %g\n",diffen,orben,ratio);

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

double dotprod(vector a, vector b)
{
  double prod=0.0;

  prod=a.x*b.x+a.y*b.y+a.z*b.z;
  return prod;
}

double magnitude(vector a)
{
   double mag=0.0;

   mag=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
   return mag;
}
