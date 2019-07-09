/* Converts Cartesian coordinates in bt 
format to Orbital elements. Follows a 
single object. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 39.47841760435743 /* Newton's constant for AU,yr,Msun */
#define M 1.0e5    /* Mass of SMBH in solar masses */
#define PI 3.141592653589793  /* pi to unnecessary precision */


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
vector vecdiv(vector a, double d);  /* Division: a/d. */
vector crossprod(vector a, vector b); /* Cross product, axb */
double magnitude(vector a);  /* Magnitude of vector a */

main()
{

  int Obj,N_steps,Out_int,nDigits;
   printf("Enter object number  to follow: ");
   scanf("%d",&Obj);
   printf("Enter number of steps: ");
   scanf("%d",&N_steps);
   printf("Enter output digits: ");
   scanf("%d",&nDigits);
   printf("Enter output interval: ");
   scanf("%d",&Out_int);

   int i,j,k;
   double j0,j1,j2;
   double L,minternal;
   vector *r,*v;
   double j10,j11,j12,j13,energy;
   double a,e,incl;
   vector angmom,amomhat;
   char fp[50],fout[50],buf[5000];
   FILE *output,*pkd;

   sprintf(fout,"1elem%d.dat",Obj);
   printf("%s\n",fout);
   output=fopen(fout,"w");
   zerovec(&angmom);

   sprintf (fp,"ss.%0*d.bt",nDigits,Out_int);
   pkd=fopen(fp,"r");
   if (!fp){
     printf("File cannot be opened\n");
     return EXIT_FAILURE;
   }
   int N_obj = 0;
   while((fgets(buf,5000,pkd)) != NULL) {
     N_obj++;
   }
   printf("%d stars\n",N_obj);
   double m[N_obj],dist[N_obj];
   r=calloc(N_obj,sizeof(vector));
   v=calloc(N_obj,sizeof(vector));
   fclose(pkd);

  for (i=Out_int; i<=N_steps; i+=Out_int)
    {
      sprintf (fp,"ss.%0*d.bt",nDigits,i);
      printf("%s\n",fp);
      pkd=fopen(fp,"r");

      for (j=0; j<N_obj; j++)
	{
	  fscanf(pkd,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
		 &j0,&j1,&m[j],&j2,&r[j].x,&r[j].y,&r[j].z,&v[j].x,&v[j].y,&v[j].z,&j10,&j11,&j12,&j13);
	  dist[j] = magnitude(r[j]);
	}
      fclose(pkd);

      minternal = M;
      for (j=0; j<N_obj; j++)
	{
	  if (dist[Obj] > dist[j]) {
	    minternal+=m[j];
	  }
	}
      zerovec(&angmom);
      angmom=vecadd(angmom,vecmult(m[Obj],crossprod(r[Obj],v[Obj])));
      amomhat = vecdiv(angmom,magnitude(angmom));
      incl=acos(amomhat.z)*180.0/PI;
      L=magnitude(angmom)*2*PI;
      energy=0.5*m[Obj]*magnitude(v[Obj])*magnitude(v[Obj])*4*PI*PI;
      energy-=G*minternal*m[Obj]/dist[Obj];
      a=-G*minternal*m[Obj]/(2.0*energy);
      e=1.0-L*L/(m[Obj]*m[Obj]*minternal*G*a);
      e=sqrt(e);
      fprintf(output,"%d %lg %lg %lg %lg %lg\n",i,a,e,incl,a*(1-e),dist[Obj]);
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
   vector sum;

   sum.x=a.x-b.x;
   sum.y=a.y-b.y;
   sum.z=a.z-b.z;
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
