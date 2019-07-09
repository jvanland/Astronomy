/* Takes NBODY6 initial conditions and
   recenters the center of mass and momentum */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
   double x;     /* x component of vector. */
   double y;     /* y component of vector. */
   double z;     /* z component of vector. */
} vector;

void zerovec(vector *r);    /* Set components of vector to zero */
vector vecadd(vector a, vector b);   /* Vector addition a+b. */
vector vecsub(vector a, vector b);   /* Vector subtraction a-b. */
vector vecmult(double d, vector a);  /* Multiplication: da. */
vector vecdiv(double d, vector a);  /* Division: a/d. */
double magnitude(vector a);  /* Magnitude of vector a */

main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Infile\n", argv[0]);
    return 1;
  }

  int i;
  double mtot;
  vector *r,*v;
  vector mom,com;
  char *Infile = argv[1];
  char fout[50],buf[5000];
  FILE *input,*output;

  input=fopen(Infile,"r");
  sprintf(fout,"%s.center",Infile);
  output=fopen(fout,"w");

  int N_obj = 0;
  while((fgets(buf,5000,input)) != NULL) {
    N_obj++;
  }
  rewind(input);
  printf("N_obj = %d\n",N_obj);
  double m[N_obj];

  r=calloc(N_obj,sizeof(vector));
  v=calloc(N_obj,sizeof(vector));
  zerovec(&mom);
  zerovec(&com);
  mtot=0;

  for (i=0; i<N_obj; i++)
    {
	  fscanf(input,"%lg %lg %lg %lg %lg %lg %lg",
		 &m[i],&r[i].x,&r[i].y,&r[i].z,&v[i].x,&v[i].y,&v[i].z);
	  mom=vecadd(mom,vecmult(m[i],v[i]));
	  com=vecadd(com,vecmult(m[i],r[i]));
	  mtot+=m[i];
    }

  mom = vecdiv(mtot,mom);
  com = vecdiv(mtot,com);
  printf("com offset: %f mom offset: %f\n",magnitude(com),magnitude(mom));

  for (i=0; i<N_obj; i++)
    {
      r[i] = vecsub(r[i],com);
      v[i] = vecsub(v[i],mom);
      fprintf(output,"%lg %lg %lg %lg %lg %lg %lg\n",
	      m[i],r[i].x,r[i].y,r[i].z,v[i].x,v[i].y,v[i].z);
    }

  printf("wrote %s\n",fout);

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
   vector sub;

   sub.x=a.x-b.x;
   sub.y=a.y-b.y;
   sub.z=a.z-b.z;
   return sub;
}

vector vecmult(double d, vector a)
{
   vector mult;

   mult.x=d*a.x;
   mult.y=d*a.y;
   mult.z=d*a.z;
   return mult;
}

vector vecdiv(double d, vector a)
{
   vector div;

   div.x=a.x/d;
   div.y=a.y/d;
   div.z=a.z/d;
   return div;
}

double magnitude(vector a)
{
   double mag=0.0;

   mag=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
   return mag;
}
