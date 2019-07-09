#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

int main()
{

  int i,j[100];

  for (i = 0; i < 100; i++)
    {
      if (i > 10) {
	j[i] = i; }
      else j[i] = 10*i;
    }

  
  
}
    
