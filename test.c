#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

int main()
{
  int N_init = 1956,lost = 4;
  char cmd[100];

  sprintf(cmd,"awk 'FNR == 2 {$1 = %d}1' intry > tmp && mv tmp intry",N_init-lost);
  system(cmd);

}
