/* Restart nbody6 every time it stalls out, until
   the run is complete. Sends an email warning. */

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define NAME "10kruns_1" /* Name of version */
#define TIME 2507.0 /* Total runtime */
#define TSTART 214.62 /* Runtime completed */
#define INTERVAL 2 /* Output steps to retreat from stall */
#define CNT2 10 /* Restart iteration */

int main()
{

  int sz1,sz2,cnt1,cnt2=7,lost,IPID,N_init,N_out,N_res;
  float runtime,freq,T_end1,T_end2,T_start,T_tot,DTMIN,RMIN;
  char cmd[100],buf[5000];
  FILE *output,*intry,*vars,*fp,*mail;

  // Keep going until run is complete 
  T_tot = TSTART;
  cnt2 = CNT2;
  T_end1 = 100;
  system("./nbody6 < intry > output &"); 
  system("echo $! > kill.dat");
  fp = fopen("kill.dat","r");
  fscanf(fp,"%d",&IPID);
  fclose(fp);
  while(T_tot < TIME) {
    mail = fopen("mail.txt","w");
    fprintf(mail,"Project %s restarted!\n",NAME);
    cnt2++;
    fprintf(mail,"cnt2 = %d\n",cnt2);
    /* Check every 10 min. that "output" has changed size,
     else exit loop  and kill nbody6 */
    sz1 = 3;
    sz2 = 5;
    cnt1 = 0;
    while(sz1 != sz2) {
      //Manual exit condition
      if (fp = fopen("EXIT","r"))  {
        output = fopen("exit.dat","w");
	fprintf(output,"nbodydriver in %s exited on manual intervention\n",NAME);
	fclose(output);
	fclose(fp);
        system("rm -f EXIT");
	return 1;
      }
      cnt1++;
      output = fopen("output","rb");
      fseek(output, 0L, SEEK_END);
      sz2 = ftell(output);
      fclose(output);
      sz1 = sz2;
      system("sleep 10m");
      output = fopen("output","rb");
      fseek(output, 0L, SEEK_END);
      sz2 = ftell(output);
      fclose(output);
    }
    sprintf(cmd,"kill %d",IPID);
    system(cmd);
    fprintf(mail,"Stalled size = %d, runtime = %d minutes\n",sz2,cnt1*10);

  /* Create new initial conditions file slightly
     before nbody6 stalled */
    system("awk 'FNR == 2 {print $1}' intry > vars.dat");
    system("awk 'FNR == 3 {print $4}' intry >> vars.dat");
    system("awk 'FNR == 3 {print $6}' intry >> vars.dat");
    system("tac output | awk '$1==\"ADJUST:\"{print $4;exit}' >> vars.dat");
    system("awk 'FNR == 9 {print $1}' intry >> vars.dat");
    system("awk 'FNR == 9 {print $2}' intry >> vars.dat");
    vars = fopen("vars.dat","r");
    fscanf(vars,"%d",&N_init);
    fscanf(vars,"%f",&freq);
    fscanf(vars,"%f",&runtime);
    fscanf(vars,"%f",&T_end2);
    fscanf(vars,"%f",&DTMIN);
    fscanf(vars,"%f",&RMIN);
    fclose(vars);
    N_out = (int) (T_end2/freq);
    N_res = N_out - INTERVAL;
    T_start = T_end2 - INTERVAL*freq;
    if(T_start <= 0.0) T_start = 0.0;
    fprintf(mail,"N_init = %d, N_out = %d, N_res = %d\n",N_init,N_out,N_res);
    sprintf(cmd,"./restart_fort %d %d %d %d",N_init,N_out,N_res,N_init);
    system(cmd);

    // Move old files to new
    if(T_start > 0.0 || T_end1 < 1.0) {
      sprintf(cmd,"mv OUT3 OUT3_%d",cnt2);
      system(cmd);
      sprintf(cmd,"mv fort.10 fort.10_%d",cnt2);
      system(cmd);
      sprintf(cmd,"mv output output_%d",cnt2);
      system(cmd);
      sprintf(cmd,"cp intry intry_%d",cnt2);
      system(cmd);
      system("mv fort.restart fort.10");
      // Count # of escapers
      lost = 0;
      if (fp = fopen("ESC","r")) {
	while((fgets(buf,5000,fp)) != NULL) {
	  lost++;
	}
	fclose(fp);
	sprintf(cmd,"mv ESC ESC_%d",cnt2);
	system(cmd);
      }
      fprintf(mail,"N_esc = %d\n",lost);
    }
    else {
      cnt2 --;
      system("rm -f OUT3 output fort.restart");
      if (fp = fopen("ESC","r")) {
    	system("rm -f ESC");
    	fclose(fp);
      }
      fprintf(mail,"cnt2 = %d\n",cnt2);
    }
    // Remove unnecessary files
    if (fp = fopen("fort.1","r")) { 
    	system("rm -f fort.1"); 
    	fclose(fp);
    }
    if (fp = fopen("OUT9","r")) { 
    	system("rm -f OUT9"); 
    	fclose(fp);
    } 
    system("rm -f fort.2 OUT33");

    //Exit Conditions
    T_tot += T_start;
    if(T_tot + 1.0 > TIME) {
      fp = fopen("exit.dat","w");
      fprintf(fp,"nbodydriver exited in %s because the run is complete\n",NAME);
      fclose(fp);
      system("mailx -s \"Exit\" jvanland@astro.umd.edu < exit.dat");
      return 1;
    }
    if(T_end1 <= 1.0 && T_end2 <= 1.0) {
      fp = fopen("exit.dat","w");
      fprintf(fp,"nbodydriver exited in %s because restart won't help\n",NAME);
      fclose(fp);
      system("mailx -s \"Exit\" jvanland@astro.umd.edu < exit.dat");
      return 1;
    }

    // Edit "intry"
    fprintf(mail,"T_end2 = %f, T_tot = %f, TIME = %f\n",T_end2,T_tot,TIME);
    sprintf(cmd,"awk 'FNR == 2 {$1 = %d}1' intry > tmp && mv tmp intry",N_init-lost);
    system(cmd);
    sprintf(cmd,"awk 'FNR == 3 {$6 = %f}1' intry > tmp && mv tmp intry",TIME-T_tot);
    system(cmd);
    // Change regularization parameters if not running smoothly
    if(T_end2 <= 20.0) {
      if(DTMIN == 0.0) {
	system("awk 'FNR == 9 {$1 = 1.0E-11}1' intry > tmp && mv tmp intry");
	system("awk 'FNR == 9 {$2 = 1.0E-3}1' intry > tmp && mv tmp intry");
      }
      else {
	system("awk 'FNR == 9 {$1 = 0.0E-11}1' intry > tmp && mv tmp intry");
	system("awk 'FNR == 9 {$2 = 0.0E-3}1' intry > tmp && mv tmp intry");
      }
      fprintf(mail,"DTMIN & RMIN adjusted\n");
    }
    fclose(mail);

    //Send email with details
    system("mailx -s \"Restart\" jvanland@astro.umd.edu < mail.txt");
    system ("cp mail.txt info.dat");

    //Restart nbody6
    system("./nbody6 < intry > output &");
    T_end1 = T_end2;
    system("echo $! > kill.dat");
    fp = fopen("kill.dat","r");
    fscanf(fp,"%d",&IPID);
    fclose(fp);

  }

}
