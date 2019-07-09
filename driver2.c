/* Driver code to follow progress of binary along already
established path of its center of mass.  Will the Kozai
effect play any role? 
Starts from fulllist_check.dat
Requires:
pkdgrav output file
hndrag (inputPNwen.hnb, output.hnb, drag.in)
reposition.c
Output is fulllist.dat (process with process2.c)
*/

#include <stdio.h>
#include <stdlib.h>


main(int argc, char *argv[ ])
{
  if(argc != 2) {
    printf("Usage: %s Data_file\n", argv[0]);
    return 1;
  }
   char *pkdata = argv[1];
   char buf[5000],cmd[100];
   int i,t,N=0,status;
   double xsuper,ysuper,zsuper,vxsuper,vysuper,vzsuper;

   FILE *pkd,*super;

   pkd=fopen(pkdata,"r");

   while((fgets(buf,5000,pkd)) != NULL) {
     N++;
   }
   rewind(pkd);

   fscanf(pkd,"%d %lg %lg %lg %lg %lg %lg",
      &t,&xsuper,&ysuper,&zsuper,&vxsuper,&vysuper,&vzsuper);
   fscanf(pkd,"%d %lg %lg %lg %lg %lg %lg",
      &t,&xsuper,&ysuper,&zsuper,&vxsuper,&vysuper,&vzsuper);
   super=fopen("super.dat","w");
   fprintf(super,"%d %.16g %.16g %.16g %.16g %.16g %.16g",
      t,xsuper,ysuper,zsuper,vxsuper,vysuper,vzsuper);
   fclose(super);
   system("head -2 fulllist_check.dat > fulllist.dat");
   system("head -2 fulllist_check.dat > initial.dat");
   system("./reposition");
   system("hndrag inputPNwen.hnb");

   for (i=1; i<N; i++)
   {
      fscanf(pkd,"%d %lg %lg %lg %lg %lg %lg",
         &t,&xsuper,&ysuper,&zsuper,&vxsuper,&vysuper,&vzsuper);
      super=fopen("super.dat","w");
      fprintf(super,"%d %.16g %.16g %.16g %.16g %.16g %.16g",
         t,xsuper,ysuper,zsuper,vxsuper,vysuper,vzsuper);
      fclose(super);
      system("tail -n1 body1.dat >> fulllist.dat");
      system("tail -n1 body2.dat >> fulllist.dat");
      system("tail -n1 body1.dat > initial.dat");
      system("tail -n1 body2.dat >> initial.dat");
      system("./reposition");
      system("hndrag inputPNwen.hnb");
   }
   fclose(pkd);
   system("tail -n1 body1.dat >> fulllist.dat");
   system("tail -n1 body2.dat >> fulllist.dat");
   sprintf(cmd,"echo %s complete! >> ~/LUSTRE/nbody6/runs/done.dat",pkdata);
   system(cmd);
}
