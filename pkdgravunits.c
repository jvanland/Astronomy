/* Converts data from nbody6 format, back to pkdgrav
format, with the use of convert.dat */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793
#define Rconv 0.000440724
#define Vconv 1.0601
#define Mtot 2019

int main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Initfile\n", argv[0]);
    return 1;
  }

   int N_obj = 0;
   char buf[5000],fout[50];
   char *Initfile = argv[1];
   FILE *initial,*output;
   sprintf(fout,"%s.convert",Initfile);
   output=fopen(fout,"w");
   initial=fopen(Initfile,"r");
   while((fgets(buf,5000,initial)) != NULL) {
     N_obj++;
   }
   rewind(initial);
   printf("N_obj = %d\n",N_obj);

   int i,j;
   double m,x,y,z;
   double vx,vy,vz;

   /* Load in data and convert*/
   for (i=0; i<N_obj; i++)
   {
      fscanf(initial,"%lg %lg %lg %lg %lg %lg %lg",
         &m,&x,&y,&z,&vx,&vy,&vz);
      m *= Mtot;
      x /= Rconv;
      y /= Rconv;
      z /= Rconv;
      vx /= Vconv;
      vy /= Vconv;
      vz /= Vconv;
      //  fprintf(output,"%d %d %lg 0.01 %lg %lg %lg %lg %lg %lg 0 0 0 3\n",i,i,m,x,y,z,vx,vy,vz);
      fprintf(output,"%lg %lg %lg %lg %lg %lg %lg\n",m,x,y,z,vx,vy,vz);
   }
   fclose(initial);

}
