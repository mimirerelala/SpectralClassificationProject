//Determines the spectral type based on the Ca II K line
#include <stdio.h>
#include <math.h>
#include "util.h"
void sp2class(float spcode, float lumcode, float metcode, 
	      float *x, float *y, int *n);

float sptCaK(float spcode)
{
  extern float *X,*Y;
  extern float Lumcode;
  extern int N;
  float *x,*y;
  float chi;
  float ssum,psum;
  float chi2 = 0.0;
  int i,n;

  if(spcode < 13.0) return(10.0);

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,Lumcode,0.0,x,y,&n);

  ssum = psum = 0.0;

  for(i=0;i<n;i++) {
    if(x[i] >= 3918.0 && x[i] <= 3925.0) {
      ssum += y[i];
      psum += Y[i];
    }
  }
  ssum /= 8.0;
  psum /= 8.0;

  if(spcode < 13.0) chi2 = 10.0;
  else {
    for(i=0;i<n;i++) {
      if(x[i] >= 3927.0 && x[i] <= 3937.0) {
        chi = y[i]/ssum - Y[i]/psum;
        chi2 += chi*chi;
      }
    }
  }
  free_vector(x,0,N);
  free_vector(y,0,N);
  return(chi2);
} 
