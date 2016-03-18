// Calculates a chi2 difference between the program star and a given library spectrum
#include <stdio.h>
#include <math.h>
#include "util.h"
void sp2class(float spcode, float lumcode, float metcode, 
	      float *x, float *y, int *n);

float spt2min(float pcode[])
{
  extern float *X,*Y;
  extern int N;
  extern float wlow;
  extern float whigh;
  extern double scool;
  extern double shot;
  float *x,*y;
  float chi;
  float chi2 = 0.0;
  int i,n;
  float spcode,lumcode;

  spcode = pcode[1];
  lumcode = pcode[2];

  if(spcode >= scool) return(1.0e+5);
  if(spcode <= shot) return(1.0e+5);

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);

  for(i=0;i<n;i++) {
    if(x[i] >= wlow+100.0 && x[i] <= whigh-100.0) {
      chi = y[i] - Y[i];
      chi2 += chi*chi;
    }
  }
  chi2 /= (whigh-wlow-200.0);

  free_vector(x,0,N);
  free_vector(y,0,N);
  return(chi2);
}
