/* Converts data from pkdgrav units to nbody units,
both of which have G = 1 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793
#define RUNTIME 112000.0 /* Runtime in yrs */

int main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Initfile\n", argv[0]);
    return 1;
  }

   int N_obj = 0;
   char buf[5000],fout[50];
   char *Initfile = argv[1];
   FILE *initial,*output,*convert;
   //   sprintf(fout,"%s.convert",Initfile);
   output=fopen("fort.10","w");
   convert=fopen("convert.dat","w");
   initial=fopen(Initfile,"r");
   while((fgets(buf,5000,initial)) != NULL) {
     N_obj++;
   }
   rewind(initial);
   fprintf(convert,"N_obj = %d\n",N_obj);

   int i,j,var,runtime;
   double m[N_obj],x[N_obj],y[N_obj],z[N_obj],d;
   double vx[N_obj],vy[N_obj],vz[N_obj];
   double kinetic, potential,mtot,qvir,q,beta;

   /* Load in data and calculate Q */
   kinetic = 0.0;
   for (i=0; i<N_obj; i++)
   {
      fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg",
         &m[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i]);
      mtot +=m[i];
     kinetic+=0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
   }
   fclose(initial);
   fprintf(convert,"Mtot = %lg\n",mtot);
   potential=0.0;
   for (i=0; i<N_obj-1; i++)
   {
      for (j=i+1; j<N_obj; j++)
      {
         d=sqrt((x[j]-x[i])*(x[j]-x[i])+(y[j]-y[i])*(y[j]-y[i])+(z[j]-z[i])*(z[j]-z[i]));
         potential-=m[i]*m[j]/d;
      }
   }
   q = -kinetic/potential;
   if(q < 0.45 || q > 0.55) var = 1;
   else var = 2;
   fprintf(convert,"kinetic = %f, potential = %f, Q = %f\n",kinetic, potential, q);
   fprintf(convert,"%d\n",var);
   /* Rescale masses and velocities and recalculate energies */
   kinetic=0.0;
   for (i=0; i<N_obj; i++)
   {
     m[i] = m[i]/mtot;
     vx[i] = vx[i]/sqrt(mtot);
     vy[i] = vy[i]/sqrt(mtot);
     vz[i] = vz[i]/sqrt(mtot);
     kinetic+=0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
   }
   potential=0.0;
   for (i=0; i<N_obj-1; i++)
   {
      for (j=i+1; j<N_obj; j++)
      {
         d=sqrt((x[j]-x[i])*(x[j]-x[i])+(y[j]-y[i])*(y[j]-y[i])+(z[j]-z[i])*(z[j]-z[i]));
         potential-=m[i]*m[j]/d;
      }
   }
   /* Rescale positions and velocities with Beta */
   beta = -(potential+kinetic)/0.25;
   qvir = 1/sqrt(mtot);
   fprintf(convert,"Beta = %lg, Qvir = %lg\n",beta,qvir);
   fprintf(convert,"Rconv = %lg, Vconv = %lg, Tconv = %lg\n",beta,qvir/sqrt(beta),beta*sqrt(beta)/qvir);
   kinetic = 0.0;
   for (i=0; i<N_obj; i++)
     {
       x[i]*=beta;
       y[i]*=beta;
       z[i]*=beta;
       vx[i] = vx[i]/sqrt(beta);
       vy[i] = vy[i]/sqrt(beta);
       vz[i] = vz[i]/sqrt(beta);
       kinetic+=0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
       fprintf(output,"%lg %lg %lg %lg %lg %lg %lg\n",
	      m[i],x[i],y[i],z[i],vx[i],vy[i],vz[i]);
     }
   potential = 0.0;
   for (i=0; i<N_obj-1; i++)
   {
      for (j=i+1; j<N_obj; j++)
      {
         d=sqrt((x[j]-x[i])*(x[j]-x[i])+(y[j]-y[i])*(y[j]-y[i])+(z[j]-z[i])*(z[j]-z[i]));
         potential-=m[i]*m[j]/d;
      }
   }
   runtime = (int) RUNTIME*2.0*PI*beta*sqrt(beta)/qvir;
   fprintf(convert,"Etot = %lg, kinetic = %lg, potential = %lg\n",kinetic+potential,kinetic,potential);
   fprintf(convert,"runtime = %d 1 pc = %lg\n",runtime,1.0/(2.063e5*beta));
}
