/* Initial condition generator.  We assume that a supermassive
black hole of mass M solar masses has N stars around it of
mass m.  These stars have a number density n~r^{-alpha}, e.g.,
alpha=7/4 gives the Bahcall-Wolf distribution.  We assume that
the stars are distributed in a spherically symmetric way, and
that the velocity distribution is isotropic at all points, with
a magnitude drawn from a Maxwell-Boltzmann distribution with
a mean equal to the velocity of a circular orbit at that
location, including the SMBH mass and the masses of the stars
interior to the orbit (this guarantees that the virial theorem
holds).  Output is in units of AU for distance, and AU/yr for
the three velocity components. 

In this version, we prevent the total velocity squared from
being more than 1.9 or less than 0.1 times the circular
velocity squared at a given location.  This maintains energy symmetry,
hence allows us to continue to obey the virial theorem, but
prevents us from having unbound stars in our initial conditions. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 1.0     /* Newton's constant for AU,yr,yr/2pi */
#define M 1.0e3   /* Mass of SMBH in solar masses */
#define N 1000      /* Number of stars */
#define Nbin 1   /* Number of binaries */
#define m 1.0      /* Mass of stars in solar masses; assumed constant */
#define mbin 20.0    /* Mass of binaries */
#define alpha 2.0  /* Slope of number density power law */
#define rmin 0.0   /* Minimum allowed radius in AU */
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
vector vecdiv(vector a, double d);  /* Division: a/d. */
vector crossprod(vector a, vector b); /* Cross product, axb */
double magnitude(vector a);  /* Magnitude of vector a */
double gasdev(void);

int main()
{
   int i,j,iPid,cnt1,cnt2;
   double rvar;
   double d[N],rinfl,rmax,ctheta,stheta,phi,mint;
   double mass,Lrel,L,a,e,energy;
   double v2,v3d2,vel;
   vector r,v,angmom;
   FILE *initial;
   
   srand48(iPid=getpid());
   initial=fopen("fort.11","w");

   rinfl = 2.063e5*sqrt(M/1.0e6);
   rmax = rinfl*(N/M);

   fprintf(initial,"%lg 0.0 0.0 0.0 0.0 0.0 0.0\n", M);

   for (i=0; i<N; i++) {
     rvar=drand48();
     d[i]=pow((pow(rmin,3.0-alpha)+rvar*(pow(rmax,3.0-alpha)-pow(rmin,3.0-alpha))),1.0/(3.0-alpha));
     if(i==N-1) d[i] = 5000.0;
   }

   for (i=0; i<N; i++)
   {
     //   do {
      ctheta=1.0-2.0*drand48();
      stheta=sqrt(1.0-ctheta*ctheta);
      phi=2.0*PI*drand48();
      r.x=d[i]*stheta*cos(phi);
      r.y=d[i]*stheta*sin(phi);
      r.z=d[i]*ctheta;
      if(i==N-1) {
	r.x = 5000.0;
	r.y = 0.0;
	r.z = 0.0;
      }
      mint=M;
       for (j=N-Nbin-1; j<N; j++) 
         {
	   if (d[i] > d[j]) {
	     mint+=mbin; //add mass of binaries if interior
           }       
         }
       for (j=0; j<N-Nbin; j++)
	 {
	   if (d[i] > d[j]) {
	     mint+=m; //add mass of stars
           }
	 }
      v3d2=G*mint/d[i];  /* 3-D vel disp squared */
      do
      {
         v2=v3d2+0.3*v3d2*gasdev();
      } while (v2<0.1*v3d2 || v2>1.9*v3d2);
      vel=sqrt(v2);
      ctheta=1.0-2.0*drand48();
      stheta=sqrt(1.0-ctheta*ctheta);
      phi=2.0*PI*drand48();
      v.x=vel*stheta*cos(phi);
      v.y=vel*stheta*sin(phi);
      v.z=vel*ctheta;
      if(i==N-1) {
	v.x = 0.0;
	v.y = sqrt(G*mint/d[i]);
	v.z = 0.0;
      }
      /* Determine a and e */
      if(i<N-1) mass = m;
      else mass = mbin;
      zerovec(&angmom);
      angmom=vecadd(angmom,vecmult(mass,crossprod(r,v)));
      L=magnitude(angmom);
      energy=0.5*mass*magnitude(v)*magnitude(v);
      energy-=G*mint*mass/d[i];
      a=-G*mint*mass/(2.0*energy);
      e=1.0-L*L/(mass*mass*mint*G*a);
      e=sqrt(e);
      Lrel = (1.0-e*e)*sqrt(a/rmax)*(mbin/m);
      //   } while (Lrel < 1.0);
     printf("%d %lg\n",i,Lrel);
      if (i<N-Nbin) {
	fprintf(initial,"%lg %lg %lg %lg %lg %lg %lg\n",m,r.x,r.y,r.z,v.x,v.y,v.z);
      }
      else {
	fprintf(initial,"%lg %lg %lg %lg %lg %lg %lg\n",mbin,r.x,r.y,r.z,v.x,v.y,v.z);
      }
   }
   fclose(initial);
}

double gasdev(void)
{
/* Generate gaussian random deviate with mean zero and variance 1. */

   static double gset;
   double fac,r,v1,v2;

     do
     {
       v1=2.0*drand48()-1.0;
       v2=2.0*drand48()-1.0;
       r=v1*v1+v2*v2;
     } while (r>=1.0 || r==0.0);
     fac=sqrt(-2.0*log(r)/r);
     gset=v1*fac;
     return v2*fac;
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
