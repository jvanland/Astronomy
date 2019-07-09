#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1000

typedef struct
{
   double x;     /* x component of vector. */
   double y;     /* y component of vector. */
   double z;     /* z component of vector. */
} vector;

 void zerovec(vector *r);    /* Set components of vector to zero */
 vector vecsub(vector a, vector b);   /* Vector subtraction a-b. */
 double magnitude(vector a);  /* Magnitude of vector a */

main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Initfile.bt\n", argv[0]);
    return 1;
  }

   char *Initfile = argv[1];
   int i;
   double id,j1,j2,j3;
   vector *r,*v;
   double j10,j11,j12,j13;
   double distance1,distance2;
   FILE *initial;

   initial=fopen(Initfile,"r");
   r=calloc(N,sizeof(vector));
   v=calloc(N,sizeof(vector));

   for (i=0; i<N; i++)
   {
     fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
       &id,&j1,&j2,&j3,&r[i].x,&r[i].y,&r[i].z,&v[i].x,&v[i].y,&v[i].z,&j10,&j11,&j12,&j13);
   }

   for (i=0; i<N-2; i++)
   {
     distance1=magnitude(vecsub(r[i],r[998]));
     distance2=magnitude(vecsub(r[i],r[999]));
     printf("%d %lg %lg\n",i,distance1,distance2);
   }

   fclose(initial);
   free(r);
   free(v);
}

void zerovec(vector *r)
{
   (*r).x=0.0;
   (*r).y=0.0;
   (*r).z=0.0;
}

vector vecsub(vector a, vector b)
{
   vector sub;

   sub.x=a.x-b.x;
   sub.y=a.y-b.y;
   sub.z=a.z-b.z;
   return sub;
}

double magnitude(vector a)
{
   double mag=0.0;

   mag=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
   return mag;
}
