#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793

main(int argc, char *argv[ ])
{
  if(argc != 3) {
    printf("Usage: %s Input Output\n", argv[0]);
    return 1;
  }
   char *Input = argv[1];
   char *Output = argv[2];
   int i,N=0;
   float h1,h2,h4,h11,h12,h13,h14;
   double m,x,y,z,vx,vy,vz;
   FILE *input,*output;
   char buf[5000];

   input=fopen(Input,"r");
   output=fopen(Output,"w");
   while((fgets(buf,5000,input)) != NULL) {
     N++;
   }
   rewind(input);

   for (i=0; i<N; i++)
   {
     fscanf(input,"%f %f %lg %f %lg %lg %lg %lg %lg %lg %f %f %f %f",
	    &h1,&h2,&m,&h4,&x,&y,&z,&vx,&vy,&vz,&h11,&h12,&h13,&h14);

     fprintf(output,"%f %.16g %.16g %.16g %.16g %.16g %.16g\n",
	     m,x,y,z,2*PI*vx,2*PI*vy,2*PI*vz);
   }


}
