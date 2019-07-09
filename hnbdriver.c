/* Driver code to follow progress of binary along already
established path of its center of mass.  Will the Kozai
effect play any role? 
Requires:
c.o.m. orbit file
hndrag (input.hnb, user.in)
binmaker.c
binplacer.c
binreplacer.c
Output is fullbody1,2.dat (process with process3.c)
*/

#include <stdio.h>
#include <stdlib.h>


main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Data_file\n", argv[0]);
    return 1;
  }
   char *comdata = argv[1];
   char buf[5000];
   int i,t,N=0;
   double xsuper,ysuper,zsuper,vxsuper,vysuper,vzsuper;

   FILE *com,*super;

   com=fopen(comdata,"r");

   while((fgets(buf,5000,com)) != NULL) {
     N++;
   }
   rewind(com);

   system("hndrag input.hnb");

   fscanf(com,"%d %lg %lg %lg %lg %lg %lg",
      &t,&xsuper,&ysuper,&zsuper,&vxsuper,&vysuper,&vzsuper);

   for (i=1; i<N; i++)
   {
      fscanf(com,"%d %lg %lg %lg %lg %lg %lg",
         &t,&xsuper,&ysuper,&zsuper,&vxsuper,&vysuper,&vzsuper);
      super=fopen("comtemp.dat","w");
      fprintf(super,"%d %.16g %.16g %.16g %.16g %.16g %.16g",
         t,xsuper,ysuper,zsuper,vxsuper,vysuper,vzsuper);
      fclose(super);
      system("tail -n1 body1.dat >> fullbody1.dat");
      system("tail -n1 body2.dat >> fullbody2.dat");
      system("tail -n1 body1.dat > bintemp.dat");
      system("tail -n1 body2.dat >> bintemp.dat");
      system("./binreplacer");
      system("hndrag input.hnb");
   }
   fclose(com);
   system("tail -n1 body1.dat >> fullbody1.dat");
   system("tail -n1 body2.dat >> fullbody2.dat");
}
